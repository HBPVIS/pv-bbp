/*=========================================================================

  Program: DevSample_PickPixelValue

  Copyright (c) Mark Wyszomierski
  All rights reserved.
  See copyright.txt or http://www.devsample.org/copyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the above copyright notice for more information.

=========================================================================*/
#include "VtkObserverErrorWarning.h"

using namespace std;

VtkObserverErrorWarning::VtkObserverErrorWarning()
{
}

void VtkObserverErrorWarning::Execute(vtkObject *caller, 
                                      unsigned long ul, 
                                      void* pData) 
{
    const char* pszDetails = (const char*)pData;
    switch (ul) {
        case vtkCommand::ErrorEvent:
            m_listErrors.push_back(pszDetails);
            break;
        case vtkCommand::WarningEvent:
            m_listWarnings.push_back(pszDetails);
            break;
    }
}

string VtkObserverErrorWarning::CreateDetailsStringFromList(const list<string>& theList) const
{
    string str;
    for (list<string>::const_iterator it  = theList.begin();
                                      it != theList.end();
                                      it++)
    {
        str += *it + "\n";
    }
    return str;
}

string VtkObserverErrorWarning::CreateDetailsString(eMsgType msgType) const
{
    string str;
    switch (msgType) {
        case ERROR:
            return CreateDetailsStringFromList(m_listErrors);
        case WARNING:
            return CreateDetailsStringFromList(m_listWarnings);
        default:
            return CreateDetailsStringFromList(m_listErrors) +
                   CreateDetailsStringFromList(m_listWarnings);
    }
}

void VtkObserverErrorWarning::Reset()
{
    m_listErrors.clear();
    m_listWarnings.clear();
}