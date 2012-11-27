#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkCommand.h"
#include "vtkPointPicker.h"
#include "vtkSphereSource.h"
#include "vtkProperty.h"

#include <stdio.h>

class vtkMyCallback : public vtkCommand
{
public:
        vtkMyCallback() { m_pvtkActorSelection = NULL; };
        static vtkMyCallback *New() { return new vtkMyCallback; }
        void PrintSelf(ostream&, vtkIndent) { }
        void PrintTrailer(ostream&, vtkIndent) { }
        void PrintHeader(ostream&, vtkIndent) { }
        void CollectRevisions(ostream&) {}
        void SetSelectionActor(vtkActor* pvtkActorSelection) { m_pvtkActorSelection = pvtkActorSelection; };

        virtual void Execute(vtkObject *caller, unsigned long, void*)
        {
                vtkRenderWindowInteractor *iren = reinterpret_cast<vtkRenderWindowInteractor*>(caller);
                vtkPointPicker *picker = (vtkPointPicker *)iren->GetPicker();
                cout << "PointId: " << picker->GetPointId() << "\n";
                if (picker->GetPointId() != -1)
                {
                        if (m_pvtkActorSelection)
                                m_pvtkActorSelection->SetPosition(picker->GetPickPosition());
                        iren->Render();
                }
        }
private:
        vtkActor* m_pvtkActorSelection;
};

int main( int argc, char *argv[] )
{
        int i;
        static float x[8][3]={{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
                {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}};
        static vtkIdType pts[6][4]={{0,1,2,3}, {4,5,6,7}, {0,1,5,4},
                {1,2,6,5}, {2,3,7,6}, {3,0,4,7}};

        // We'll create the building blocks of polydata including data attributes.
        vtkPolyData *cube = vtkPolyData::New();
        vtkPoints *points = vtkPoints::New();
        vtkCellArray *polys = vtkCellArray::New();
        vtkFloatArray *scalars = vtkFloatArray::New();

        // Load the point, cell, and data attributes.
        for (i=0; i<8; i++) points->InsertPoint(i,x[i]);
        for (i=0; i<6; i++) polys->InsertNextCell(4,pts[i]);
        for (i=0; i<8; i++) scalars->InsertTuple1(i,i);

        // We now assign the pieces to the vtkPolyData.
        cube->SetPoints(points);
        cube->SetPolys(polys);
        cube->GetPointData()->SetScalars(scalars);

        // Now we'll look at it.
        vtkPolyDataMapper *cubeMapper = vtkPolyDataMapper::New();
        cubeMapper->SetInput(cube);
        cubeMapper->SetScalarRange(0,7);
        vtkActor *cubeActor = vtkActor::New();
        cubeActor->SetMapper(cubeMapper);

        // The usual rendering stuff.
        vtkCamera *camera = vtkCamera::New();
        camera->SetPosition(1,1,1);
        camera->SetFocalPoint(0,0,0);

        vtkRenderer *renderer = vtkRenderer::New();
        vtkRenderWindow *renWin = vtkRenderWindow::New();
        renWin->AddRenderer(renderer);

        vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
        iren->SetRenderWindow(renWin);

        renderer->AddActor(cubeActor);
        renderer->SetActiveCamera(camera);
        renderer->ResetCamera();
        renderer->SetBackground(1,1,1);

        renWin->SetSize(300,300);

        // create marker for pick
        vtkSphereSource *sphere = vtkSphereSource::New();
        sphere->SetThetaResolution(8); sphere->SetPhiResolution(8);
        sphere->SetRadius(0.1);
        vtkPolyDataMapper *sphereMapper = vtkPolyDataMapper::New();
        sphereMapper->SetInput(sphere->GetOutput());
        vtkActor *sphereActor = vtkActor::New();
        sphereActor->SetMapper(sphereMapper);
        sphereActor->GetProperty()->SetColor(1,0,0);
        sphereActor->PickableOff();

        renderer->AddActor(sphereActor);

        // start renderer
        renWin->Render();

        // init picker
        vtkPointPicker *picker = vtkPointPicker::New();
        picker->SetTolerance(0.01);
        iren->SetPicker(picker);

        // init callback
        vtkMyCallback *callback = vtkMyCallback::New();
        callback->SetSelectionActor(sphereActor);
        iren->AddObserver(vtkCommand::EndPickEvent, callback);

        // start interaction
        iren->Start();

        iren->RemoveObserver(callback);

        // Clean up
        cube->Delete();
        points->Delete();
        polys->Delete();
        scalars->Delete();
        cubeMapper->Delete();
        cubeActor->Delete();
        camera->Delete();
        sphere->Delete();
        sphereMapper->Delete();
        sphereActor->Delete();
        picker->Delete();
        renderer->Delete();
        renWin->Delete();
        iren->Delete();
        callback->Delete();

        return 0;
}

