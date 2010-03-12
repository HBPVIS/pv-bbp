/*=========================================================================

  Program: DevSample_PickPixelValue

  Copyright (c) Mark Wyszomierski
  All rights reserved.
  See copyright.txt or http://www.devsample.org/copyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the above copyright notice for more information.

=========================================================================*/
#include "WrapperVtkDisplayDICOMImage.h"
#include "VtkObserverErrorWarning.h"
#include "VtkObserverMouseMove.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkDICOMImageReader.h"
#include "vtkCamera.h"
#include "vtkImageMapToWindowLevelColors.h"
#include "vtkImageData.h"
#include "vtkImageActor.h"

using namespace std;

WrapperVtkDisplayDICOMImage::WrapperVtkDisplayDICOMImage(int cxWin, int cyWin)
{
    m_pIren = vtkRenderWindowInteractor::New();
    m_pImageViewer = vtkImageViewer2::New();
    m_pPicker = vtkPointPicker::New();
    m_pTextInfo = WrapperVtkText::New();
    
    // Setup the basic rendering pipeline.
    SetupSceneBasics();

    ResetSize(cxWin, cyWin); // Set the size right away.
}

WrapperVtkDisplayDICOMImage::~WrapperVtkDisplayDICOMImage()
{
    m_pIren->Delete();
    m_pImageViewer->Delete();
    m_pPicker->Delete();
    m_pTextInfo->Delete();
}

void WrapperVtkDisplayDICOMImage::SetupSceneBasics()
{
    // Just attach our interactor to the image viewer.
    m_pImageViewer->SetupInteractor(m_pIren);

    // Turn interpolation of the image off so it's easier to see
    // pixel boundaries for this example.
    m_pImageViewer->GetImageActor()->InterpolateOff();

    // Connect our picker with our interactor, set the picking tolerance to zero.
    m_pPicker->SetTolerance(0.0);
    m_pIren->SetPicker(m_pPicker);

    // Add a text overlay so we can show pixel values of the mouse over.
    m_pTextInfo->AddToScene(m_pImageViewer->GetRenderer());
    m_pTextInfo->SetText("");
    m_pTextInfo->SetPosition(WrapperVtkText::POS_LOWER_LEFT);
    m_pTextInfo->Show(1);
}

void WrapperVtkDisplayDICOMImage::ResetSize(int cx, int cy)
{
    // Reset the render window size, reposition the text actors.
    m_pImageViewer->GetRenderWindow()->SetSize(cx, cy);
    m_pTextInfo->SetParentWindowSize(cx, cy);
}

bool WrapperVtkDisplayDICOMImage::LoadVolumeFromFolder(string strPath)
{
    // We need a specific handler to load different file types.
    vtkDICOMImageReader* pRead = vtkDICOMImageReader::New();
    if (!LoadVolumeFromFolder_Dicom(strPath, pRead)) {
        pRead->Delete();
        return false;
    }

    // Set the input of the image viewer to be whatever we read from the
    // DICOM directory. This won't work with true color DICOM images, they'll
    // just come out greyscale too.
    m_pImageViewer->SetInput(pRead->GetOutput());
    m_pImageViewer->SetColorWindow(255); // Set an initial window that makes sense for you.
    m_pImageViewer->SetColorLevel(128);  // Set an initial level  that makes sense for you.
    m_pImageViewer->GetRenderer()->ResetCamera();

    // Add a mouse move obvserver which we will have tell us the pixel value
    // under the mouse whenever it moves.
    VtkObserverMouseMove *observeMouseMove = VtkObserverMouseMove::New(m_pImageViewer, m_pIren, m_pPicker, m_pTextInfo->GetTextMapper());
    m_pIren->AddObserver(vtkCommand::MouseMoveEvent, observeMouseMove);
    observeMouseMove->Delete();

    // Get rid of the DICOM directory reader.
    pRead->Delete();

    return true;
}

bool WrapperVtkDisplayDICOMImage::LoadVolumeFromFolder_Dicom(string strPath, 
                                                        vtkDICOMImageReader* pRead)
{
    // Read all the DICOM files in the specified directory.
    pRead->SetDirectoryName(strPath.c_str());

    // Watch for errors and warnings while loading the images.
    VtkObserverErrorWarning *observeError = VtkObserverErrorWarning::New();
    pRead->AddObserver(vtkCommand::ErrorEvent,   observeError);
    pRead->AddObserver(vtkCommand::WarningEvent, observeError);

    // Ok start reading all the files in the directory now.
    pRead->Update();
    if (observeError->DidErrorOccur()) {
        m_strDetails = "Error reading DICOM images:\n" + observeError->CreateDetailsString(VtkObserverErrorWarning::ERROR);
        observeError->Delete();
        return false;
    }
    observeError->Delete();

    return true;
}

void WrapperVtkDisplayDICOMImage::Start()
{
    // win32 applications should override this call.
    // Initialize() makes the render window appear, and if done
    // before SetParentId(HWND) will make the render window appear
    // outside the CWnd frame! Start() should not be called at all
    // for win32 applications, it begins the render window running
    // in a new separate thread.
    m_pIren->Initialize();
	m_pIren->Start();
}