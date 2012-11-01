/*=========================================================================

   Program: ParaView
   Module:    $RCSfile: pqCircuitReaderPanel.cxx,v $

   Copyright (c) 2005-2008 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaView is a free software; you can redistribute it and/or modify it
   under the terms of the ParaView license version 1.2. 

   See License_v1.2.txt for the full ParaView license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "pqCircuitReaderPanel.h"
#include "ui_pqCircuitReaderPanel.h"

// Qt includes
#include <QAction>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QTimer>
#include <QTreeWidget>
#include <QVariant>
#include <QVector>
#include <QMap>
#include <QHeaderView>

// VTK includes

// ParaView Server Manager includes
#include "vtkEventQtSlotConnect.h"
#include "vtkGraph.h"
#include "vtkProcessModule.h"
#include "vtkPVArrayInformation.h"
#include "vtkPVDataInformation.h"
#include "vtkPVDataSetAttributesInformation.h"
#include "vtkPVSILInformation.h"
#include "vtkSmartPointer.h"
#include "vtkSMProperty.h"
#include "vtkSMProxyManager.h"
#include "vtkSMSourceProxy.h"
#include "vtkSMPropertyHelper.h"

// ParaView includes
#include "pqPropertyManager.h"
#include "pqProxy.h"
#include "pqServer.h"
#include "pqTimeKeeper.h"
#include "pqSMAdaptor.h"
#include "pqTreeWidget.h"
#include "pqTreeWidgetSelectionHelper.h"
#include "pqTreeWidgetItemObject.h"
#include "vtkSMDoubleVectorProperty.h"
#include "pqSILModel.h"
#include "pqProxySILModel.h"
#include "pqTreeViewSelectionHelper.h"
#include "pqTreeView.h"
//----------------------------------------------------------------------------
class pqCircuitReaderPanel::pqUI : public QObject, public Ui::vtkCircuitReaderPanel 
{
public:
  pqUI(pqCircuitReaderPanel* p) : QObject(p) 
   {
   this->VTKConnect = vtkSmartPointer<vtkEventQtSlotConnect>::New();
   this->SILUpdateStamp = -1;
   }

  pqSILModel SILModel;
  QVector<double> TimestepValues;
  QMap<QTreeWidgetItem*, QString> TreeItemToPropMap;
  vtkSmartPointer<vtkEventQtSlotConnect> VTKConnect;
  int SILUpdateStamp;
};
//----------------------------------------------------------------------------
pqCircuitReaderPanel::pqCircuitReaderPanel(pqProxy* object_proxy, QWidget* p) :
  Superclass(object_proxy, p)
{
  this->UI = new pqUI(this);

  // Create a frame and copy our custom gui into it
  QFrame *frame = new QFrame();
  frame->setObjectName(QString::fromUtf8("frame"));
  frame->setFrameShape(QFrame::StyledPanel);
  frame->setFrameShadow(QFrame::Raised);
  frame->setFrameStyle(QFrame::NoFrame);
  this->UI->setupUi(frame);
  // add custom panel to the existing auto-generated stuff
  int rows = this->PanelLayout->rowCount();
  int cols = this->PanelLayout->columnCount();
  QVBoxLayout* subLayout = new QVBoxLayout();
  subLayout->addWidget(frame, 1);
  subLayout->setMargin(0);
  subLayout->setSpacing(0);

  this->PanelLayout->addLayout(subLayout, rows-1, 0, -1, -1);
  //
  this->DisplItem = 0;
  //
  this->UI->XMLFileName->setServer(this->referenceProxy()->getServer());
  //
  this->UI->VTKConnect->Connect(
    this->referenceProxy()->getProxy(),
    vtkCommand::UpdateInformationEvent,
    this, SLOT(updateSIL()));

  pqProxySILModel* blockProxyModel = new pqProxySILModel("DataSet Types", &this->UI->SILModel);
  blockProxyModel->setSourceModel(&this->UI->SILModel);
  this->UI->Blocks->setModel(blockProxyModel);
  this->UI->Blocks->header()->setClickable(true);
  QObject::connect(this->UI->Blocks->header(), SIGNAL(sectionClicked(int)),
    blockProxyModel, SLOT(toggleRootCheckState()), Qt::QueuedConnection);
  //
  pqProxySILModel* dimsProxyModel = new pqProxySILModel("Dummy Name", &this->UI->SILModel);
  dimsProxyModel->setSourceModel(&this->UI->SILModel);
  QString title("Dimensions of DataSet");
  dimsProxyModel->setHeaderTitle(title);
  this->UI->Dimensions->setModel(dimsProxyModel);
  this->UI->Dimensions->header()->setClickable(false);
  this->UI->Dimensions->header()->setStretchLastSection(true);

  this->updateSIL();
  //
  this->UI->Blocks->expandAll();
  this->UI->Blocks->header()->setStretchLastSection(true);
  blockProxyModel->setData(QModelIndex(), Qt::Unchecked, Qt::CheckStateRole);
  
  this->linkServerManagerProperties();


//  QObject::connect(this->UI->Blocks,
//    SIGNAL(clicked(const QModelIndex &)),
//    this, SLOT(blockItemChanged(const QModelIndex &)));

  QObject::connect(this->UI->Blocks->selectionModel(),
    SIGNAL(currentChanged(const QModelIndex &, const QModelIndex &)),
    this, SLOT(blockItemChanged(const QModelIndex &, const QModelIndex &)));

  QList<pqTreeWidget*> treeWidgets = this->findChildren<pqTreeWidget*>();
  foreach (pqTreeWidget* tree, treeWidgets)
    {
    new pqTreeWidgetSelectionHelper(tree);
    }

  QList<pqTreeView*> treeViews = this->findChildren<pqTreeView*>();
  foreach (pqTreeView* tree, treeViews)
    {
    new pqTreeViewSelectionHelper(tree);
    }
}
//----------------------------------------------------------------------------
pqCircuitReaderPanel::~pqCircuitReaderPanel()
{
}
//----------------------------------------------------------------------------
void pqCircuitReaderPanel::updateSIL()
{
/*
  vtkSMProxy* reader = this->referenceProxy()->getProxy();
  reader->UpdatePropertyInformation(reader->GetProperty("SILUpdateStamp"));

  int stamp = vtkSMPropertyHelper(reader, "SILUpdateStamp").GetAsInt();
  if (stamp != this->UI->SILUpdateStamp)
    {
    this->UI->SILUpdateStamp = stamp;
    vtkProcessModule* pm = vtkProcessModule::GetProcessModule();
    vtkPVSILInformation* info = vtkPVSILInformation::New();
    pm->GatherInformation(reader->GetConnectionID(),
      vtkProcessModule::DATA_SERVER, info,
      reader->GetID());
    this->UI->SILModel.update(info->GetSIL());
    info->Delete();
    }
*/
}
//----------------------------------------------------------------------------
void pqCircuitReaderPanel::addSelectionsToTreeWidget(const QString& prop, 
                                      QTreeWidget* tree)
{
  vtkSMProperty* SMProperty = this->proxy()->GetProperty(prop.toAscii().data());
  QList<QVariant> SMPropertyDomain;
  SMPropertyDomain = pqSMAdaptor::getSelectionPropertyDomain(SMProperty);
  int j;
  for(j=0; j<SMPropertyDomain.size(); j++)
    {
    QString varName = SMPropertyDomain[j].toString();
    this->addSelectionToTreeWidget(varName, varName, tree, prop, j);
    }
}
//----------------------------------------------------------------------------
void pqCircuitReaderPanel::addSelectionToTreeWidget(const QString& name,
                                               const QString& realName,
                                               QTreeWidget* tree,
                                               const QString& prop,
                                               int propIdx)
{
  vtkSMProperty* SMProperty = this->proxy()->GetProperty(prop.toAscii().data());

  if(!SMProperty || !tree)
    {
    return;
    }

  QList<QString> strs;
  strs.append(name);
  pqTreeWidgetItemObject* item;
  item = new pqTreeWidgetItemObject(tree, strs);
  item->setData(0, Qt::ToolTipRole, name);
  this->propertyManager()->registerLink(item, 
                      "checked", 
                      SIGNAL(checkedStateChanged(bool)),
                      this->proxy(), SMProperty, propIdx);

  this->UI->TreeItemToPropMap[item] = prop;
}
//----------------------------------------------------------------------------
void pqCircuitReaderPanel::linkServerManagerProperties()
{
  this->propertyManager()->registerLink(
    this->UI->Blocks->model(), "values", SIGNAL(valuesChanged()),
    this->proxy(),
    this->proxy()->GetProperty("VariablesArrayStatus"));

  // parent class hooks up some of our widgets in the ui
  this->Superclass::linkServerManagerProperties();

  QObject::connect(this->UI->Refresh,
    SIGNAL(pressed()), this, SLOT(onRefresh()));
}
//----------------------------------------------------------------------------
QString pqCircuitReaderPanel::formatDataFor(vtkPVArrayInformation* ai)
{
  QString info;
  if(ai)
    {
    int numComponents = ai->GetNumberOfComponents();
    int dataType = ai->GetDataType();
    double range[2];
    for(int i=0; i<numComponents; i++)
      {
      ai->GetComponentRange(i, range);
      QString s;
      if(dataType != VTK_VOID && dataType != VTK_FLOAT && 
         dataType != VTK_DOUBLE)
        {
        // display as integers (capable of 64 bit ids)
        qlonglong min = qRound64(range[0]);
        qlonglong max = qRound64(range[1]);
        s = QString("%1 - %2").arg(min).arg(max);
        }
      else
        {
        // display as reals
        double min = range[0];
        double max = range[1];
        s = QString("%1 - %2").arg(min,0,'f',6).arg(max,0,'f',6);
        }
      if(i > 0)
        {
        info += ", ";
        }
      info += s;
      }
    }
  else
    {
    info = "Unavailable";
    }
  return info;
}


void pqCircuitReaderPanel::onRefresh()
{
  vtkSMSourceProxy* sp = vtkSMSourceProxy::SafeDownCast(this->proxy());
  vtkSMProperty *prop = sp->GetProperty("Refresh");

  // The "Refresh" property has no values, so force an update this way
  prop->SetImmediateUpdate(1);
  prop->Modified();

  // "Pull" the values
  sp->UpdatePropertyInformation(sp->GetProperty("TimeRange"));
  sp->UpdatePropertyInformation(sp->GetProperty("TimestepValues")); 
}
//----------------------------------------------------------------------------
void pqCircuitReaderPanel::blockItemChanged(const QModelIndex &current, const QModelIndex &previous)
{
  // a child leaf is a dataset, but we want the node name to fetch dimension information
  // so go up one if the user selected a leaf
  QModelIndex index = current;
  if (!this->UI->SILModel.hasChildren(index)) {
    index = this->UI->SILModel.parent(index);
  }
  QVariant value = this->UI->SILModel.data(index , Qt::DisplayRole);
  QString sval = value.toString();
  //
  pqProxySILModel* dimsProxyModel = new pqProxySILModel(sval, &this->UI->SILModel);
  dimsProxyModel->setNoCheckBoxes(true);
  dimsProxyModel->setSourceModel(&this->UI->SILModel);
  sval = "Dimensions : " + value.toString();
  dimsProxyModel->setHeaderTitle(sval);
  this->UI->Dimensions->setModel(dimsProxyModel);
  this->UI->Dimensions->header()->setClickable(false);
  this->UI->Dimensions->header()->setStretchLastSection(true);
}
//----------------------------------------------------------------------------
