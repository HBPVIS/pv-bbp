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

int main(int argc, char* argv[])
{
    // There should be one argument:
    // 1) path to folder of DICOM images
    if (argc < 2) {
        cout << "Not enough parameters:" << endl;
        cout << "Example:" << endl;
        cout << "    > DevSample_PickPixelValue.exe C:\\test" << endl << endl;
        return -1;
    }

    std::string strFolderPath = argv[1];

    // Create the wrapper and load the DICOM files.
    WrapperVtkDisplayDICOMImage w(640, 480);
    if (w.LoadVolumeFromFolder(strFolderPath)) {
        w.Start();
    }
    else {
        cout << "Error:" << endl;
        cout << w.GetDetails() << endl << endl;
    }

    return 0;
}