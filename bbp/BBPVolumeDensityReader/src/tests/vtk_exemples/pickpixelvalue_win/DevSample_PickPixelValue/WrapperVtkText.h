/*=========================================================================

  Program: DevSample_PickPixelValue

  Copyright (c) Mark Wyszomierski
  All rights reserved.
  See copyright.txt or http://www.devsample.org/copyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the above copyright notice for more information.

=========================================================================*/
#ifndef __WrapperVtkText_h
#define __WrapperVtkText_h

#include "vtkObject.h"

class vtkRenderer;
class vtkTextActor;
class vtkTextMapper;
class vtkTextProperty;
class vtkActor2D;

class WrapperVtkText : public vtkObject {

    public:
        static WrapperVtkText *New();

        // Just pass the current window size and the text will center itself.
        void SetParentWindowSize(int cxWin, int cyWin);

        // Add the to the scene.
        void AddToScene(vtkRenderer *pParentRenderer);

        // Hide or show the text.
        void Show(int bShow);

        // Set the text.
        void SetText(const char *pszText);

        // Set the position of the text.
        void SetPosition(int nPosition);

        // Get a pointer to the internal text mapper.
        vtkTextMapper* GetTextMapper() { return m_TextMapper; }

        enum {
            POS_CENTER,
            POS_UPPER_LEFT,
            POS_UPPER_RIGHT,
            POS_LOWER_RIGHT,
            POS_LOWER_LEFT
        };

    private:

		// Combine a text mapper with an actor.
		vtkTextMapper *m_TextMapper;
		vtkActor2D    *m_Actor;

        int m_cxWin, m_cyWin;
        int m_nPositionType;

    private:
        WrapperVtkText(const WrapperVtkText&); // Not implemented
        void operator=(const WrapperVtkText&); // Not implemented

    protected:
        WrapperVtkText();
        ~WrapperVtkText();
};
#endif