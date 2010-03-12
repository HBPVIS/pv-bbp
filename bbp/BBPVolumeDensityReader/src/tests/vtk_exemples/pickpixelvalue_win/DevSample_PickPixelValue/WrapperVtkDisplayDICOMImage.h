/*=========================================================================

  Program: DevSample_PickPixelValue

  Copyright (c) Mark Wyszomierski
  All rights reserved.
  See copyright.txt or http://www.devsample.org/copyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the above copyright notice for more information.

=========================================================================*/
#ifndef __WrapperVtkDisplayDICOMImage_h
#define __WrapperVtkDisplayDICOMImage_h

#include "vtkRenderWindowInteractor.h"
#include "vtkDICOMImageReader.h"
#include "vtkImageViewer2.h"
#include "vtkPointPicker.h"
#include "WrapperVtkText.h"
#include <string>
#include <sstream>

#define FmtStr(str_out, strm) { std::ostringstream oss; oss << strm; str_out = oss.str(); }

class WrapperVtkDisplayDICOMImage {

    public:
        // Create the empty render wrapper with the passed dimensions.
        WrapperVtkDisplayDICOMImage(int cxWin = 480, int cyWin = 480);
        virtual ~WrapperVtkDisplayDICOMImage();

        // Load all images in the passed folder path for volume rendering.
        bool LoadVolumeFromFolder(std::string strPath);

        // Start the interactor (only used with non-GUI apps!).
        void Start();

        // Set the size of the rendering window.
        void ResetSize(int cx, int cy);

        // Get error details etc.
        std::string GetDetails() { return m_strDetails; }

    protected:

        // Setup the basic rendering pipeline.
        void SetupSceneBasics();

        // Handle specifics for loading DICOM files from a folder.
        bool LoadVolumeFromFolder_Dicom(std::string strPath, 
                                        vtkDICOMImageReader* pReader);

    protected:

        vtkRenderWindowInteractor*         m_pIren;
        vtkImageViewer2*                   m_pImageViewer;
        vtkPointPicker*                    m_pPicker;
        WrapperVtkText*                    m_pTextInfo;
  
        // Used to record errors etc.
        std::string m_strDetails;
};
#endif