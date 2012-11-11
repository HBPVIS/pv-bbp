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

  pqProxySILModel* targetsProxyModel = new pqProxySILModel("Targets", &this->UI->SILModel);
  targetsProxyModel->setSourceModel(&this->UI->SILModel);
  this->UI->Targets->setModel(targetsProxyModel);
//  this->UI->Targets->header()->setClickable(true);
//  QObject::connect(this->UI->Targets->header(), SIGNAL(sectionClicked(int)),
//    targetsProxyModel, SLOT(toggleRootCheckState()), Qt::QueuedConnection);
  //
  this->updateSIL();
  //
  this->UI->Targets->expandAll();
  this->UI->Targets->header()->setStretchLastSection(true);
  targetsProxyModel->setData(QModelIndex(), Qt::Unchecked, Qt::CheckStateRole);
  
  this->linkServerManagerProperties();

//  QObject::connect(this->UI->Targets,
//    SIGNAL(clicked(const QModelIndex &)),
//    this, SLOT(targetItemChanged(const QModelIndex &)));

  QObject::connect(this->UI->Targets->selectionModel(),
    SIGNAL(currentChanged(const QModelIndex &, const QModelIndex &)),
    this, SLOT(targetItemChanged(const QModelIndex &, const QModelIndex &)));

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
  vtkSMProxy* reader = this->referenceProxy()->getProxy();
  reader->UpdatePropertyInformation(reader->GetProperty("SILUpdateStamp"));

  int stamp = vtkSMPropertyHelper(reader, "SILUpdateStamp").GetAsInt();
  if (stamp != this->UI->SILUpdateStamp) {
    this->UI->SILUpdateStamp = stamp;
    vtkPVSILInformation* info = vtkPVSILInformation::New();
    reader->GatherInformation(info);
    if (info->GetSIL()) {
      this->UI->SILModel.update(info->GetSIL());
    }
    info->Delete();
  }
}

//----------------------------------------------------------------------------
void pqCircuitReaderPanel::addSelectionsToTreeWidget(const QString& prop, QTreeWidget* tree)
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
    this->UI->Targets->model(), "values", SIGNAL(valuesChanged()),
    this->proxy(),
    this->proxy()->GetProperty("TargetsStatus"));

  // parent class hooks up some of our widgets in the ui
  this->Superclass::linkServerManagerProperties();
  // on startup, we want to unselect everything
  pqProxySILModel* targetsProxyModel = (pqProxySILModel*)(this->UI->Targets->model());
  targetsProxyModel->setData(QModelIndex(), Qt::Unchecked, Qt::CheckStateRole);

  // and select the default target
  vtkSMProxy* reader = this->referenceProxy()->getProxy();
  std::string default = vtkSMPropertyHelper(reader, "DefaultTarget").GetAsString();

  pqSILModel* silModel = qobject_cast<pqSILModel*>(targetsProxyModel->sourceModel());
  if (silModel->findVertex(default.c_str())!=-1) {
    targetsProxyModel->setData(silModel->makeIndex(silModel->findVertex(default.c_str())), Qt::Checked, Qt::CheckStateRole);
  }
}
//----------------------------------------------------------------------------
void pqCircuitReaderPanel::targetItemChanged(const QModelIndex &current, const QModelIndex &previous)
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
}
//----------------------------------------------------------------------------
