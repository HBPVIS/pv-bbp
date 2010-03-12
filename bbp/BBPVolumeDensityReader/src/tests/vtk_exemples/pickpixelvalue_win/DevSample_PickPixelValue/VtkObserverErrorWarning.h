/*=========================================================================

  Program: DevSample_PickPixelValue

  Copyright (c) Mark Wyszomierski
  All rights reserved.
  See copyright.txt or http://www.devsample.org/copyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the above copyright notice for more information.

=========================================================================*/
#ifndef __WrapperVtkObservers_h
#define __WrapperVtkObservers_h

#include "vtkCommand.h"
#include <list>
#include <string>

// Some VTK classes don't have an error return method, rather you need to observe
// errors. This class can be used to observe errors and warnings on an object. Any
// errors or warnings will be stored in the member lists for review.
class VtkObserverErrorWarning : public vtkCommand {

    public:
        VtkObserverErrorWarning();
        static VtkObserverErrorWarning* New() { return new VtkObserverErrorWarning(); }

        // Message type.
        enum eMsgType {
            ERROR = 0,
            WARNING,
            ALL
        };

        // Combine all recorded errors and warnings into one formatted string.
        std::string CreateDetailsString(eMsgType msgType) const;

        // Check if an error or warning occurred.
        bool DidErrorOccur()   const { return !m_listErrors.empty(); }
        bool DidWarningOccur() const { return !m_listWarnings.empty(); }

        // Clear the internal lists.
        void Reset();

    protected:

        std::list<std::string> m_listErrors;
        std::list<std::string> m_listWarnings;
        std::string CreateDetailsStringFromList(const std::list<std::string>& theList) const;

    private:

        void Execute(vtkObject *caller, unsigned long ul, void* pData);
};
#endif