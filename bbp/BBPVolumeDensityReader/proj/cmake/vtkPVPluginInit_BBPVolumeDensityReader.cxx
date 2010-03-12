// This file is generated.  Do not edit.


/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile: vtkPVPluginInit.cxx.in,v $

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#if 1
#include "/home/lasserre/workspace/paraview/dev/bbp/BBPVolumeDensityReader/proj/cmake/vtkSMXML_BBPVolumeDensityReader.h"
#endif // HAVE_XML

#ifdef _WIN32
// __cdecl gives an unmangled name
#define C_DECL __cdecl
#define C_EXPORT extern "C" __declspec(dllexport)
#else
#define C_DECL
#define C_EXPORT extern "C"
#endif

// entry point to get Plugin name a string
C_EXPORT const char* C_DECL ParaViewPluginName()
{
  return "BBPVolumeDensityReader";
}

// entry point to get Plugin version as a string
C_EXPORT const char* C_DECL ParaViewPluginVersion()
{
  return "1.0";
}

// entry point to get PluginRequiredOnServer as an int 
C_EXPORT int ParaViewPluginRequiredOnServer()
{
  return 1;
}

// entry point to get PluginRequiredOnClient as an int 
C_EXPORT int ParaViewPluginRequiredOnClient()
{
  return 1;
}

#if 0
// entry point to get Plugin-Depended-Plugins as a string
C_EXPORT const char* ParaViewPluginRequiredPlugins()
{
  return "";
}
#endif

#if 1

namespace {
  class StaticInitXML
  {
    public:
    StaticInitXML()
    {
      static char* xmls[] = 
      {
             BBPVolumeDensityReaderVolumeDensityReaderGetInterfaces()
      };

      XMLString = xmls;
      NumberOfStrings = sizeof(xmls) / sizeof(char*);
    }
    ~StaticInitXML()
    {
      // clean up new'd arrays
      for(int i=0; i<NumberOfStrings; i++)
      {
        delete [] XMLString[i];
      }       
    }
    char** XMLString;
    int NumberOfStrings;
  };

}

// entry point to get XML as a string
C_EXPORT void C_DECL ParaViewPluginXMLList(int& num, char** & xml)
{
  static StaticInitXML staticinit;
  num = staticinit.NumberOfStrings;
  xml = staticinit.XMLString;
}
#endif // HAVE_XML

#if 1

#include "vtkProcessModule.h"
#include "vtkClientServerInterpreter.h"

extern "C" void BBPVolumeDensityReader_Initialize(vtkClientServerInterpreter *arlu);

namespace {
  class StaticInitSMWrappings
    {
  public:
    StaticInitSMWrappings(void (*initfunc)(vtkClientServerInterpreter*))
      {
      // This will call (*initfunc)(processModule->GetInterpreter()) if the
      // process module has been initialized. Otherwise, it will register the
      // callback with the process module so that when it gets intialized, the
      // callback will be called. 
      // This removes the need for calling the BBPVolumeDensityReader_Initialize function
      // explicitly if an application is linking against the plugin directly.
      vtkProcessModule::InitializeInterpreter(initfunc);
      }
    };
}
static StaticInitSMWrappings InitSMWrappings(BBPVolumeDensityReader_Initialize);

#endif // HAVE_SRCS

