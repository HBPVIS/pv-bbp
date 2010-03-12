/*=========================================================================

  Program: DevSample_PickPixelValue

  Copyright (c) Mark Wyszomierski
  All rights reserved.
  See copyright.txt or http://www.devsample.org/copyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the above copyright notice for more information.

=========================================================================*/
#include "WrapperVtkText.h"
#include "vtkRenderer.h"
#include "vtkObjectFactory.h"
#include "vtkTextActor.h"
#include "vtkTextMapper.h"
#include "vtkTextProperty.h"
#include "vtkActor2D.h"

vtkStandardNewMacro(WrapperVtkText);

WrapperVtkText::WrapperVtkText()
{
	m_TextMapper = vtkTextMapper::New();
    m_TextMapper->SetInput("");

	vtkTextProperty *pProperty = m_TextMapper->GetTextProperty();
	pProperty->SetBold(1);
    pProperty->SetItalic(0);
    pProperty->SetFontSize(14);
    pProperty->SetFontFamily(VTK_ARIAL);
    pProperty->SetJustification(VTK_TEXT_CENTERED);
    pProperty->SetVerticalJustification(VTK_TEXT_CENTERED);
    m_nPositionType = POS_CENTER;
    m_cxWin = 0;
    m_cyWin = 0;

	m_Actor = vtkActor2D::New();
	m_Actor->SetMapper(m_TextMapper);
}

WrapperVtkText::~WrapperVtkText()
{
    m_TextMapper->Delete();
	m_Actor->Delete();
}

void WrapperVtkText::AddToScene(vtkRenderer *pParentRenderer)
{
    pParentRenderer->AddActor(m_Actor);
    Show(0);
    Show(1);
}

void WrapperVtkText::SetParentWindowSize(int cxWin, int cyWin)
{
    m_cxWin = cxWin;
    m_cyWin = cyWin;
	SetPosition(m_nPositionType);
}

void WrapperVtkText::SetText(const char *pszText)
{
    m_TextMapper->SetInput(pszText);
}

void WrapperVtkText::Show(int bShow)
{
    m_TextMapper->GetTextProperty()->SetOpacity(bShow);
}

void WrapperVtkText::SetPosition(int nPosition)
{
    m_nPositionType = nPosition;
    switch (m_nPositionType) {
        case POS_CENTER:
            m_Actor->SetPosition(m_cxWin/2, m_cyWin/2);
            m_TextMapper->GetTextProperty()->SetJustification(VTK_TEXT_CENTERED);
            break;
        case POS_UPPER_LEFT:
            m_Actor->SetPosition(20, m_cyWin-20);
            m_TextMapper->GetTextProperty()->SetJustification(VTK_TEXT_LEFT);
            break;
        case POS_UPPER_RIGHT:
            m_Actor->SetPosition(m_cxWin-20, m_cyWin-20);
            m_TextMapper->GetTextProperty()->SetJustification(VTK_TEXT_RIGHT);
            break;
        case POS_LOWER_RIGHT:
            m_Actor->SetPosition(m_cxWin-20, 20);
            m_TextMapper->GetTextProperty()->SetJustification(VTK_TEXT_RIGHT);
            break;
        case POS_LOWER_LEFT:
            m_Actor->SetPosition(20, 20);
            m_TextMapper->GetTextProperty()->SetJustification(VTK_TEXT_LEFT);
            break;
    }
}