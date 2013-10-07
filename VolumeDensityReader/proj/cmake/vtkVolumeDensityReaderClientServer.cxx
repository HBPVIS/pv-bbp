// ************ AUTO GENERATED ****************
// ClientServer wrapper for vtkVolumeDensityReader object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkSystemIncludes.h"
#include "vtkVolumeDensityReader.h"
#include "vtkClientServerInterpreter.h"
#include "vtkClientServerStream.h"

#include "vtkParse.h"

//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_New(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 2)
    {
    vtkVolumeDensityReader  *temp20;
      {
      temp20 = (op)->New();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  return 0;
}
//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_GetClassName(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 2)
    {
    const char    *temp20;
      {
      temp20 = (op)->GetClassName();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  return 0;
}
//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_IsA(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 3)
    {
    char    *temp0;
    int      temp20;
    if(msg.GetArgument(0, 2, &temp0))
      {
      temp20 = (op)->IsA(temp0);
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  return 0;
}
//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_NewInstance(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 2)
    {
    vtkVolumeDensityReader  *temp20;
      {
      temp20 = (op)->NewInstance();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  return 0;
}
//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_SafeDownCast(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 3)
    {
    vtkObject  *temp0;
    vtkVolumeDensityReader  *temp20;
    if(vtkClientServerStreamGetArgumentObject(msg, 0, 2, &temp0, "vtkObject"))
      {
      temp20 = (op)->SafeDownCast(temp0);
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << (vtkObjectBase *)temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  return 0;
}
//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_GetWholeExtent(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 2)
    {
    int     *temp20;
      {
      temp20 = (op)->GetWholeExtent();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << vtkClientServerStream::InsertArray(temp20,6) << vtkClientServerStream::End;
      return 1;
      }
    }
  return 0;
}
//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_SetFileName(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 3)
    {
    char    *temp0;
    if(msg.GetArgument(0, 2, &temp0))
      {
      op->SetFileName(temp0);
      return 1;
      }
    }
  return 0;
}
//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_GetFileName(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 2)
    {
    char    *temp20;
      {
      temp20 = (op)->GetFileName();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  return 0;
}
//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_SetTarget(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 3)
    {
    char    *temp0;
    if(msg.GetArgument(0, 2, &temp0))
      {
      op->SetTarget(temp0);
      return 1;
      }
    }
  return 0;
}
//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_GetTarget(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 2)
    {
    char    *temp20;
      {
      temp20 = (op)->GetTarget();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  return 0;
}
//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_SetSectionType(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 3)
    {
    char    *temp0;
    if(msg.GetArgument(0, 2, &temp0))
      {
      op->SetSectionType(temp0);
      return 1;
      }
    }
  return 0;
}
//------------------------------------------------------------------------auto
int vtkVolumeDensityReader_GetSectionType(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream)
{
  (void)resultStream;
  if(msg.GetNumberOfArguments(0) == 2)
    {
    char    *temp20;
      {
      temp20 = (op)->GetSectionType();
      resultStream.Reset();
      resultStream << vtkClientServerStream::Reply << temp20 << vtkClientServerStream::End;
      return 1;
      }
    }
  return 0;
}

#ifndef VTK_METHOD_MAP
#include <vtkstd/map>
typedef int (*funPtr)(const vtkClientServerStream& msg, vtkVolumeDensityReader *op, vtkClientServerStream& resultStream);
typedef vtkstd::map <vtkstd::string , funPtr> vtkMethodMap;
#endif

//-------------------------------------------------------------------------auto
/*
 * vtkVolumeDensityReaderMethodMap function creates a map with key as the name of the function
 * and value as the funptr which can be directly called.
 *
 * @return the map which returns the funptr given the name
 */
vtkMethodMap& vtkVolumeDensityReaderMethodMap()
{
  static bool once = 1;
  static vtkMethodMap map;
  if(once)
    {
    once = 0;
    map["New"]=vtkVolumeDensityReader_New;
    map["GetClassName"]=vtkVolumeDensityReader_GetClassName;
    map["IsA"]=vtkVolumeDensityReader_IsA;
    map["NewInstance"]=vtkVolumeDensityReader_NewInstance;
    map["SafeDownCast"]=vtkVolumeDensityReader_SafeDownCast;
    map["GetWholeExtent"]=vtkVolumeDensityReader_GetWholeExtent;
    map["SetFileName"]=vtkVolumeDensityReader_SetFileName;
    map["GetFileName"]=vtkVolumeDensityReader_GetFileName;
    map["SetTarget"]=vtkVolumeDensityReader_SetTarget;
    map["GetTarget"]=vtkVolumeDensityReader_GetTarget;
    map["SetSectionType"]=vtkVolumeDensityReader_SetSectionType;
    map["GetSectionType"]=vtkVolumeDensityReader_GetSectionType;
    }
  return map;
}


vtkObjectBase *vtkVolumeDensityReaderClientServerNewCommand()
{
  return vtkVolumeDensityReader::New();
}

int vtkImageAlgorithmCommand(vtkClientServerInterpreter*, vtkObjectBase*, const char*, const vtkClientServerStream&, vtkClientServerStream& resultStream);

int VTK_EXPORT vtkVolumeDensityReaderCommand(vtkClientServerInterpreter *arlu, vtkObjectBase *ob, const char *method, const vtkClientServerStream& msg, vtkClientServerStream& resultStream)
{
  vtkVolumeDensityReader *op = vtkVolumeDensityReader::SafeDownCast(ob);
  if(!op)
    {
    vtkOStrStreamWrapper vtkmsg;
    vtkmsg << "Cannot cast " << ob->GetClassName() << " object to vtkVolumeDensityReader.  "
           << "This probably means the class specifies the incorrect superclass in vtkTypeRevisionMacro.";
    resultStream.Reset();
    resultStream << vtkClientServerStream::Error
                 << vtkmsg.str() << 0 << vtkClientServerStream::End;
    return 0;
    }
  (void)arlu;

  if(funPtr f = vtkVolumeDensityReaderMethodMap()[method])
  if(f(msg,op,resultStream))
    return 1;


  if (vtkImageAlgorithmCommand(arlu, op,method,msg,resultStream))
    {
    return 1;
    }
  if(resultStream.GetNumberOfMessages() > 0 &&
     resultStream.GetCommand(0) == vtkClientServerStream::Error &&
     resultStream.GetNumberOfArguments(0) > 1)
    {
    /* A superclass wrapper prepared a special message. */
    return 0;
    }
  vtkOStrStreamWrapper vtkmsg;
  vtkmsg << "Object type: vtkVolumeDensityReader, could not find requested method: \""
         << method << "\"\nor the method was called with incorrect arguments.\n";
  resultStream.Reset();
  resultStream << vtkClientServerStream::Error
               << vtkmsg.str() << vtkClientServerStream::End;
  vtkmsg.rdbuf()->freeze(0);
  return 0;
}

void vtkObject_Init(vtkClientServerInterpreter* csi);
void vtkImageAlgorithm_Init(vtkClientServerInterpreter* csi);

//-------------------------------------------------------------------------auto
void VTK_EXPORT vtkVolumeDensityReader_Init(vtkClientServerInterpreter* csi)
{
  static bool once;
  if(!once)
    {
    once = true;
    vtkObject_Init(csi);
    vtkImageAlgorithm_Init(csi);
    csi->AddNewInstanceFunction("vtkVolumeDensityReader", vtkVolumeDensityReaderClientServerNewCommand);
    csi->AddCommandFunction("vtkVolumeDensityReader", vtkVolumeDensityReaderCommand);
    }
}
