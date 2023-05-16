#include <vtkActor.h>
#include <vtkActor2D.h>  //11.12
#include <vtkCaptionActor2D.h>  //11.12
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkCellPicker.h>  //09.12
#include <vtkDataSetMapper.h>  //09.12
#include <vtkExtractSelection.h>  //09.12
#include <vtkIdTypeArray.h>  //09.12
#include <vtkInteractorStyleTrackballCamera.h>  //09.12
#include <vtkLookupTable.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>  //09.12
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkSelection.h>  //09.12
#include <vtkSelectionNode.h>  //09.12
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>  //09.12

#include <sstream>

#include <vtksys/SystemTools.hxx>
#include "vtkTimerCallback2.h"



void vtkTimerCallback2::Execute(vtkObject* caller, unsigned long eventId,
  void* vtkNotUsed(callData))
{
  vtkRenderWindowInteractor* iren = dynamic_cast<vtkRenderWindowInteractor*>(caller);
  if (vtkCommand::TimerEvent == eventId)
  {
    ++this->TimerCount;
  }
  if (TimerCount < TickLim)
  {
    k += videoTick;
    int numCells = polyData->GetNumberOfCells();
    for (int i = 0; i < numCells; i++)
    {
      int num = (int)((*colors)[k][i]) / bandPiece;
      num = num < NBand - 1 ? num : NBand - 1;
      float rgb[3] = { 255, 255 * (1 - (float)num / NBand),
                      255 * (1 - (float)num / NBand) };
      if ((*colors)[k][i] <= 0)
      {
        rgb[0] = rgb[1] = 0;
        rgb[2] = 255;
      }
      cellData->InsertTuple(i, rgb);
    }
    polyData->GetCellData()->SetScalars(cellData);
    polyData->Modified();
    mapper->Update();

    actor->RotateZ(angleA);

    iren->GetRenderWindow()->Render();
  }
  else
  {
    vtkSmartPointer<vtkMatrix4x4> M = vtkSmartPointer<vtkMatrix4x4>::New();
    M->Identity();
    actor->PokeMatrix(M);
    iren->GetRenderWindow()->Render();
    iren->DestroyTimer();
  }
}


MouseInteractorStyle::MouseInteractorStyle()
{
  selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
  selectedActor = vtkSmartPointer<vtkActor>::New();
  captionActor = vtkSmartPointer<vtkCaptionActor2D>::New();

}

void MouseInteractorStyle::OnLeftButtonDown()
{

  // Get the location of the click (in window coordinates)
  int* pos = this->GetInteractor()->GetEventPosition();
  double coords[2];

  vtkSmartPointer<vtkCellPicker> picker =
    vtkSmartPointer<vtkCellPicker>::New();
  picker->SetTolerance(0.0005);

  // Pick from this location.
  picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

  double* worldPosition = picker->GetPickPosition();
  int id = picker->GetCellId();
  if (id != -1)
  {
    //std::cout << "Cell id is: " << id << '\t' << "Value is " << colors->at(id) << std::endl;
    std::stringstream ss;
    ss << colors->at(id);
    //std::cout << "Adding number: " << ss.str() << std::endl;
    coords[0] = (double)pos[0];
    coords[1] = (double)pos[1];
    // Create an actor for the text
    captionActor->SetCaption(ss.str().c_str());
    captionActor->SetAttachmentPoint(worldPosition);
    captionActor->BorderOff();
    /* captionActor->GetCaptionTextProperty()->BoldOff();
     captionActor->GetCaptionTextProperty()->ItalicOff();
     captionActor->GetCaptionTextProperty()->ShadowOff();
    */ captionActor->ThreeDimensionalLeaderOff();

    this->CurrentRenderer->AddViewProp(captionActor);

    //std::cout << "Pick position is: " << worldPosition[0] << " "
    //    << worldPosition[1] << " " << worldPosition[2] << endl;

    vtkSmartPointer<vtkIdTypeArray> ids =
      vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);
    ids->InsertNextValue(id);

    vtkSmartPointer<vtkSelectionNode> selectionNode =
      vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::CELL);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);

    vtkSmartPointer<vtkSelection> selection =
      vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);

    vtkSmartPointer<vtkExtractSelection> extractSelection =
      vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputData(0, this->Data);
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    // In selection
    vtkSmartPointer<vtkUnstructuredGrid> selected =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extractSelection->GetOutput());

    selectedMapper->SetInputData(selected);
    selectedActor->SetMapper(selectedMapper);
    selectedActor->GetProperty()->EdgeVisibilityOn();
    //   selectedActor->GetProperty()->SetColor(
    //       colors->GetColor3d("Red").GetData());

    selectedActor->GetProperty()->SetLineWidth(10);

    this->Interactor->GetRenderWindow()
      ->GetRenderers()
      ->GetFirstRenderer()
      ->AddActor(selectedActor);
  }
  // Forward events
  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}