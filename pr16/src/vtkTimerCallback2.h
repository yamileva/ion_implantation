#pragma once

#include <vtkActor.h>
#include <vtkActor2D.h>  //11.12
#include <vtkCamera.h>
#include <vtkCaptionActor2D.h>  //11.12
#include <vtkCellData.h>
#include <vtkDataSetMapper.h>  //09.12
#include <vtkInteractorStyleTrackballCamera.h>  //09.12
#include <vtkLookupTable.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>

#include <vtksys/SystemTools.hxx>


class vtkTimerCallback2 : public vtkCommand
{
public:
  vtkTimerCallback2() = default;
  ~vtkTimerCallback2() = default;

  int timerId = 0;
  static vtkTimerCallback2* New()
  {
    vtkTimerCallback2* cb = new vtkTimerCallback2;
    cb->TimerCount = 0;
    return cb;
  }
  
  virtual void Execute(vtkObject* caller, unsigned long eventId,
                       void* vtkNotUsed(callData));
  void setVideoTick(int a) { videoTick = a; }
  void setAngleA(int g) { angleA = g; }
  void setTickLim(int al) { TickLim = al; }
  void setBandPiece(double bp) { bandPiece = bp; }
  void setNBand(int n) { NBand = n; }
  void setInitTick(int it) { k = it; }
  void setSC(int ix, int iy, int iz)  { x = ix; y = iy; z = iz; }

private:
  int TimerCount = 0;
  int videoTick = 1;
  int angleA;
  int turnV = 30;
  int TickLim;
  double bandPiece;
  int NBand;
  int k;
  int x;
  int y;
  int z;

public:
  vtkActor* actor;
  vtkPolyDataMapper* mapper;
  vtkPolyData* polyData;

  vtkUnsignedCharArray* cellData;
  std::vector<std::vector<double> >* colors;
};

class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
  static MouseInteractorStyle* New();
  MouseInteractorStyle();
  virtual void OnLeftButtonDown() override;

  vtkSmartPointer<vtkPolyData> Data;
  vtkSmartPointer<vtkDataSetMapper> selectedMapper;
  vtkSmartPointer<vtkActor> selectedActor;
  vtkSmartPointer<vtkCaptionActor2D> captionActor;

  std::vector<double>* colors;
};



