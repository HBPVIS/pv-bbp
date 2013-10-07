// Generated file.  Do not edit.

/*=========================================================================

   Program: ParaView
   Module:    $RCSfile: pqPluginImplementation.cxx.in,v $

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

========================================================================*/

#include "BBPVolumeDensityReaderPluginImplementation.h"

// we need to define QT_NO_DEBUG for release builds
// so the plugin verification data is correct for Q_EXPORT_PLUGIN2
#if defined(NDEBUG) && !defined(QT_NO_DEBUG)
#define QT_NO_DEBUG
#endif
#include <QtPlugin>



BBPVolumeDensityReaderPlugin::BBPVolumeDensityReaderPlugin()
{
  
}

BBPVolumeDensityReaderPlugin::~BBPVolumeDensityReaderPlugin()
{
}

QObjectList BBPVolumeDensityReaderPlugin::interfaces()
{
  return this->Interfaces;
}

// entry point to get Plugin name a string
const char* BBPVolumeDensityReaderPlugin::ParaViewPluginName()
{
  return "BBPVolumeDensityReader";
}

// entry point to get Plugin version as a string
const char* BBPVolumeDensityReaderPlugin::ParaViewPluginVersion() 
{
  return "1.0";
}

// entry point to get PluginRequiredOnServer as an int 
int BBPVolumeDensityReaderPlugin::ParaViewPluginRequiredOnServer()
{
  return 1;
}

// entry point to get PluginRequiredOnClient as an int 
int BBPVolumeDensityReaderPlugin::ParaViewPluginRequiredOnClient()
{
  return 1;
}

#if 0
// entry point to get Plugin-Depended-Plugins as a string
const char* BBPVolumeDensityReaderPlugin::ParaViewPluginRequiredPlugins()
{
  return "";
}
#endif

Q_EXPORT_PLUGIN2(BBPVolumeDensityReaderPlugin, BBPVolumeDensityReaderPlugin)

