/*=========================================================================

  Program: DevSample_PickPixelValue

  Copyright (c) Mark Wyszomierski
  All rights reserved.
  See copyright.txt or http://www.devsample.org/copyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the above copyright notice for more information.

=========================================================================*/
#include "VtkObserverMouseMove.h"
#include "vtkImageData.h"

using namespace std;

#define FmtStr(str_out, strm) { std::ostringstream oss; oss << strm; str_out = oss.str(); }

VtkObserverMouseMove::VtkObserverMouseMove(vtkImageViewer2* pImageViewer,
                                           vtkRenderWindowInteractor* pIren,
                                           vtkPointPicker* pPicker,
                                           vtkTextMapper* pTextMapper)
{
    m_pImageViewer = pImageViewer;
    m_pIren = pIren;
    m_pPicker = pPicker;
    m_pTextMapper = pTextMapper;
}

void VtkObserverMouseMove::Execute(vtkObject *caller, unsigned long, void*) 
{
    // Do the pick. It will return a non-zero value if we intersected the image.
    if (!m_pPicker->Pick(m_pIren->GetEventPosition()[0], 
                         m_pIren->GetEventPosition()[1], 
                         0,  // always zero.
                         m_pImageViewer->GetRenderer())) 
    {
        m_pTextMapper->SetInput("Mouse is outside of image extent.");
        m_pIren->Render();
        return;
    }

    // Get the mapped position of the mouse using the picker.
    double ptMapped[3];
    m_pPicker->GetMapperPosition(ptMapped);

    // We have to manually set the physical z-coordinate which requires
    // us to get the volume spacing.
    double dSpacing[3];
    m_pImageViewer->GetInput()->GetSpacing(dSpacing);
    ptMapped[2] = m_pImageViewer->GetSlice() * dSpacing[2];

    // Get a shortcut to the pixel data.
    vtkImageData* pImageData = m_pImageViewer->GetInput();

    // Get the volume index within the entire volume now.
    int nVolIdx = pImageData->FindPoint(ptMapped);

    // We have to handle different number of scalar components.
    switch (pImageData->GetNumberOfScalarComponents()) {
        case 1:
            {
                // Get a shortcut to the raw pixel data. This should be further
                // updated to check if your data is signed or not, but for this
                // example we'll just assume it's unsigned. You should also check
                // the type of your data too (unsigned char, unsigned short, etc).
                unsigned short* pPix = (unsigned short*)pImageData->GetScalarPointer();
                unsigned short usPix = pPix[nVolIdx];
                FmtStr(m_strDetails, "Pixel val = [" << usPix << "]   at index [" << nVolIdx << "],  coordinates(" << ptMapped[0] << "," << ptMapped[1] << ")");
            }
            break;
        
        case 3:
            {
                // For vtkImageData with multiple components, you have to get each
                // component separately. Here's one case where we assume the data
                // type is unsigned char. Add additional specific cases.
                unsigned char* pPix = (unsigned char*)pImageData->GetScalarPointer();
                unsigned char usPixR = pPix[nVolIdx * 3 + 0];
                unsigned char usPixG = pPix[nVolIdx * 3 + 1];
                unsigned char usPixB = pPix[nVolIdx * 3 + 2];
                FmtStr(m_strDetails, "Pixel val = [" << (int)usPixR << "," << (int)usPixG << "," << (int)usPixB << "]   at index [" << nVolIdx << "],  coordinates(" << ptMapped[0] << "," << ptMapped[1] << ")");
            }
            break;

        default:
            FmtStr(m_strDetails, "Unsupported number of scalar components [" << pImageData->GetNumberOfScalarComponents() << "] for pixel picker.");
            break;
    }

    // Update the info text with pixel coordinates/value if requested.
    if (m_pTextMapper) {
        m_pTextMapper->SetInput(m_strDetails.c_str());
        m_pIren->Render();
    }
}