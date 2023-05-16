#include <map>
#include <windows.h>
#include <sstream>

#include <vtkActor.h>
#include <vtkActor2D.h>  //11.12
#include <vtkCamera.h>
#include <vtkCaptionActor2D.h>  //11.12
//#include <vtkCellArray.h>  //09.12
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkCellPicker.h>  //09.12
#include <vtkDataSetMapper.h>  //09.12
#include <vtkExtractSelection.h>  //09.12
#include <vtkIdTypeArray.h>  //09.12
#include <vtkInteractorStyleTrackballCamera.h>  //09.12
#include <vtkLookupTable.h>
//#include <vtkNamedColors.h>  //09.12
//#include <vtkObjectFactory.h>  //09.12
//#include <vtkPlaneSource.h>  //09.12
//#include <vtkPoints.h>  //09.12
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>  //09.12
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkSelection.h>  //09.12
#include <vtkSelectionNode.h>  //09.12
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkUnstructuredGrid.h>  //09.12


#include <vtksys/SystemTools.hxx>

#define VIDEOMODE true

// Catch mouse events
class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
  static MouseInteractorStyle* New();

  MouseInteractorStyle()
  {
    selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    selectedActor = vtkSmartPointer<vtkActor>::New();
    captionActor = vtkSmartPointer<vtkCaptionActor2D>::New();

  }

  virtual void OnLeftButtonDown() override
  {

    // Get the location of the click (in window coordinates)
    int* pos = this->GetInteractor()->GetEventPosition();
    double coords[] = { 10., 10. };

    vtkSmartPointer<vtkCellPicker> picker =
      vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.0005);

    // Pick from this location.
    picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

    //double* worldPosition = picker->GetPickPosition();
    int id = picker->GetCellId();
    if (id != -1)
    {
      std::cout << "Cell id is: " << id << '\t' << "Value is " << colors->at(id) << std::endl;
      std::stringstream ss;
      ss << colors->at(id);
      std::cout << "Adding number: " << ss.str() << std::endl;

      // Create an actor for the text
      captionActor->SetCaption(ss.str().c_str());
      captionActor->SetAttachmentPoint(coords);
      captionActor->BorderOff();
      captionActor->GetCaptionTextProperty()->BoldOff();
      captionActor->GetCaptionTextProperty()->ItalicOff();
      captionActor->GetCaptionTextProperty()->ShadowOff();
      captionActor->ThreeDimensionalLeaderOff();

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

  vtkSmartPointer<vtkPolyData> Data;
  vtkSmartPointer<vtkDataSetMapper> selectedMapper;
  vtkSmartPointer<vtkActor> selectedActor;
  vtkSmartPointer<vtkCaptionActor2D> captionActor;

  std::vector<double>* colors;
};

vtkStandardNewMacro(MouseInteractorStyle);



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
  void setVideoTick(int a)
  {
    videoTick = a;
  }

  void setAngleA(int g)
  {
    angleA = g;
  }

  void setAngleV(int g)
  {
    angleV = g;
  }

  void setTickLim(int al)
  {
    TickLim = al;
  }

  void setBandPiece(double bp)
  {
    bandPiece = bp;
  }

  void setNBand(int n)
  {
    NBand = n;
  }

  void setInitTick(int it)
  {
    k = it;
  }

  void setSC(int ix, int iy, int iz)
  {
    x = ix; y = iy; z = iz;
  }

private:
  int TimerCount = 0;
  int videoTick = 1;
  int angleA;
  int angleV;
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
//___________________________________________________________
//--------------------------MAIN-----------------------------
//___________________________________________________________
//________________________calculate__________________________
void calculate(int angleLim, int angleTempLim, int numCells, double angleVert, double angleAxis,
  std::vector <std::vector<double> >& interCells,
  vtkTransformPolyDataFilter* transformFilter, vtkTransform* regTransform)
{
  const int rayCountX = 11;
  const int rayCountY = 41;
  const int rayCount = rayCountX * rayCountY * 2;

  double rayWidthTotal = 100, rayHeightTotal = 400;
  double rayWidth = rayWidthTotal / (rayCountX - 1),
    rayHeight = rayHeightTotal / (rayCountY - 1);

  double rayOrigin[3] = { -rayWidthTotal / 2, 150, 600 };
  double axis[3] = { 0.0, 0.0, 1.0 },
    temp[3];

  double rayStart[rayCount][3], rayEnd[rayCount][3];

  auto cellLocator = vtkSmartPointer<vtkCellLocator>::New();
  cellLocator->SetDataSet(transformFilter->GetOutput());



  std::vector <std::vector<double> > interCellsTemp(angleTempLim, std::vector<double>(numCells, 0));
  /// Данные об облученных ячейках: "доля" облучения (точная на расчетном участке)

  for (int j = 0; j < rayCountY; j++)
  {
    for (int i = 0; i < rayCountX; i++)
    {

      rayStart[i + j * rayCountX][0] = rayOrigin[0] + i * rayWidth;
      rayStart[i + j * rayCountX][1] = rayOrigin[1] + j * rayHeight;
      rayStart[i + j * rayCountX][2] = rayOrigin[2];

      rayEnd[i + j * rayCountX][0] = rayStart[i + j * rayCountX][0];
      rayEnd[i + j * rayCountX][1] = rayStart[i + j * rayCountX][1];
      rayEnd[i + j * rayCountX][2] = -rayStart[i + j * rayCountX][2];
    }
  }
  for (int j = rayCountY; j < 2 * rayCountY; j++)
  {
    for (int i = 0; i < rayCountX; i++)
    {
      rayStart[i + j * rayCountX][0] = rayOrigin[2];
      rayStart[i + j * rayCountX][1] = -rayOrigin[1] - j * rayHeight;
      rayStart[i + j * rayCountX][2] = -rayOrigin[0] - i * rayWidth;

      rayEnd[i + j * rayCountX][0] = -rayStart[i + j * rayCountX][0];
      rayEnd[i + j * rayCountX][1] = rayStart[i + j * rayCountX][1];
      rayEnd[i + j * rayCountX][2] = rayStart[i + j * rayCountX][2];
    }
  }

  for (int i = 0; i < rayCount; i++)
  {
    temp[0] = rayStart[i][0] * cos(45 / 180. * 3.141592653589793) + rayStart[i][2] * sin(45 / 180. * 3.141592653589793);
    temp[2] = -rayStart[i][0] * sin(45 / 180. * 3.141592653589793) + rayStart[i][2] * cos(45 / 180. * 3.141592653589793);
    rayStart[i][0] = temp[0];
    rayStart[i][2] = temp[2];
    temp[0] = rayEnd[i][0] * cos(45 / 180. * 3.141592653589793) + rayEnd[i][2] * sin(45 / 180. * 3.141592653589793);
    temp[2] = -rayEnd[i][0] * sin(45 / 180. * 3.141592653589793) + rayEnd[i][2] * cos(45 / 180. * 3.141592653589793);
    rayEnd[i][0] = temp[0];
    rayEnd[i][2] = temp[2];
  }


  std::map <int, std::pair<int, double> > interCellsCurrent; /// Данные об облученных ячейках: количество попаданий и накопленный cos угла атаки
  std::cout << "Идет расчет..." << std::endl;
  unsigned int start = GetTickCount();

  double xyz[3], t, pcoords[3], attack, cellNormal[3];
  double angleSaveMoment = (angleLim - 1) % angleTempLim;
  int subId, hit = 0;
  vtkIdType cellId = -1;
  int k, kTemp, kSave;
  for (k = 0, kTemp = 1, kSave = 0; k < angleLim; k++, kTemp++)  /// цикл по "времени"
  {
    cellLocator->BuildLocator();
    hit = 0;

    for (int i = 0; i < rayCount; i++)  /// цикл по точкам источника
    {
      hit = cellLocator->IntersectWithLine(rayStart[i], rayEnd[i], 0.0001, t, xyz, pcoords, subId, cellId);

      if (hit)
      {
        transformFilter->GetOutput()->GetCellData()->GetArray("Normals")->GetTuple(cellId, cellNormal);
        if (1 - std::abs(cellNormal[0] * axis[0] +
          cellNormal[1] * axis[1] +
          cellNormal[2] * axis[2]) > 1e-5)
        {
          attack = (rayEnd[i][0] - rayStart[i][0]) * cellNormal[0] +
            (rayEnd[i][1] - rayStart[i][1]) * cellNormal[1] +
            (rayEnd[i][2] - rayStart[i][2]) * cellNormal[2];
          attack /= sqrt((rayEnd[i][0] - rayStart[i][0]) * (rayEnd[i][0] - rayStart[i][0]) +
            (rayEnd[i][1] - rayStart[i][1]) * (rayEnd[i][1] - rayStart[i][1]) +
            (rayEnd[i][2] - rayStart[i][2]) * (rayEnd[i][2] - rayStart[i][2]));
          attack /= sqrt(cellNormal[0] * cellNormal[0] +
            cellNormal[1] * cellNormal[1] +
            cellNormal[2] * cellNormal[2]);

          auto it = interCellsCurrent.find(cellId);
          if (it != interCellsCurrent.end())
          {
            interCellsCurrent[cellId].first += 1;
            interCellsCurrent[cellId].second += std::abs(attack);
          }
          else
          {
            interCellsCurrent[cellId] = std::make_pair(1, std::abs(attack));
          }
        }
        else
        {
          auto it = interCellsCurrent.find(cellId);
          if (it != interCellsCurrent.end())
          {
            interCellsCurrent[cellId].first = 1;
            interCellsCurrent[cellId].second = -1;
          }
          else
          {
            interCellsCurrent[cellId] = std::make_pair(1, -1);
          }
        }
      } // конец условия наличия пересечения

      /// поворот точки источника вокруг верикальной оси
      temp[0] = rayStart[i][0] * cos(angleVert / 180. * 3.141592653589793) + rayStart[i][2] * sin(angleVert / 180. * 3.141592653589793);
      temp[2] = -rayStart[i][0] * sin(angleVert / 180. * 3.141592653589793) + rayStart[i][2] * cos(angleVert / 180. * 3.141592653589793);
      rayStart[i][0] = temp[0];
      rayStart[i][2] = temp[2];
      temp[0] = rayEnd[i][0] * cos(angleVert / 180. * 3.141592653589793) + rayEnd[i][2] * sin(angleVert / 180. * 3.141592653589793);
      temp[2] = -rayEnd[i][0] * sin(angleVert / 180. * 3.141592653589793) + rayEnd[i][2] * cos(angleVert / 180. * 3.141592653589793);
      rayEnd[i][0] = temp[0];
      rayEnd[i][2] = temp[2];
    } //конец цикла по точкам

    /// Добавление среднего cos угла атаки как доли облучения

    if (k % angleTempLim != angleSaveMoment)
    {
      for (int i = 0; i < numCells; i++)
        interCellsTemp.at(kTemp)[i] = interCellsTemp.at(kTemp - 1)[i];
      for (auto& item : interCellsCurrent)
        interCellsTemp.at(kTemp)[item.first] += item.second.second / item.second.first;
    }
    else
    {
      //      std::cout << "shift!";
      //      std::cout << k << " " << kTemp << " " << kSave << std::endl;

      for (int i = 0; i < numCells; i++)
        interCellsTemp.at(0)[i] = interCellsTemp.at(kTemp - 1)[i];
      for (auto& item : interCellsCurrent)
        interCellsTemp.at(0)[item.first] += item.second.second / item.second.first;
      for (int i = 0; i < numCells; i++)
        interCells.at(kSave)[i] = interCellsTemp.at(0)[i];

      kSave++;
      kTemp = 0;
    }
    interCellsCurrent.clear();

    if ((k + 1) % (angleLim / 10) == 0)
    {
      std::cout << (k + 1) * 100 / angleLim << "%..";
    }

    //  flyRenderWindow->Render();

      /// поворот диска вокруг своей оси
    regTransform->RotateWXYZ(angleAxis, axis);
    transformFilter->Update();



  } // конец цикла по k

//  std::cout << "end!";
//  std::cout << k << " " << kTemp << " " << kSave << std::endl;

  unsigned int end = GetTickCount();
  std::cout << std::endl << "Длительность расчета " << end - start << " ms" << std::endl;

}
//___________________________________________________________
//________________________readFile___________________________
int readFile(std::string fileName, int numCells, std::vector <std::vector<double> >& interCells)
{
  std::cout << "Чтение файла...  ";
  unsigned int start = GetTickCount();
  std::ifstream resIn(fileName);
  /* if (!resIn)
   {
     std::cout << "Ошибка чтения файла" << std::endl;
     system("pause");
     return -1;
   }*/

  int angleSaveLim, k, i;
  double num;
  resIn >> angleSaveLim >> k;
  resIn >> k >> i; resIn >> num; resIn >> num; resIn >> num;
  for (int k = 0; k < angleSaveLim; k++)
  {
    for (int i = 0; i < numCells; i++)
    {
      interCells[k][i] = 0;
    }
  }
  while (!resIn.eof())
  {
    resIn >> k >> i >> num;
    interCells[k][i] = num;
  }
  resIn.close();
  unsigned int end = GetTickCount();
  std::cout << "Чтение завершено" << std::endl;
  std::cout << "Длительность чтения " << end - start << " ms" << std::endl;
  return 0;
}

//___________________________________________________________
//________________________readInfo___________________________
int readInfo(std::string fileName, int& angleSaveLim, int& angleTempLim, int& videoAngleA, int& videoAngleV,
  double& turnVert, double& turnAxis, double& koefAxisToVert)
{


  return 0;
}

//___________________________________________________________
//________________________visualize__________________________
void visualize(double bandMin, int bandMax, double bandPiece, int Nband, double koeff, int numCells,
  std::vector <std::vector<double> >& interCells, vtkTransformPolyDataFilter* transformFilter,
  int angleSaveLim, double videoAngleA, double videoAngleV, int videoTick)
{
  auto cellData = vtkSmartPointer<vtkUnsignedCharArray>::New();
  cellData->SetNumberOfComponents(3);
  cellData->SetNumberOfTuples(numCells);

  int kInit = 0;
  if (!VIDEOMODE) {
    kInit = interCells.size() - 1;
  }
  //std::cout << "Input timestep in " << angleLim << std::endl;
  //std::cin >> kInit;
  //kInit = (angleLim - 1) % videoTick;
  //vtkSmartPointer<vtkPolyData> polyData = transformFilter->GetOutput();

  for (int i = 0; i < numCells; i++)
  {
    int num = (int)(interCells[kInit][i] - bandMin) / bandPiece;
    num = num < Nband - 1 ? num : Nband - 1;
    float rgb[3] = { 255, 255 * (1 - (float)num / Nband),
                         255 * (1 - (float)num / Nband) };
    if (interCells[kInit][i] <= 0)
    {
      rgb[0] = rgb[1] = 0;
      rgb[2] = 255;
    }
    cellData->InsertTuple(i, rgb);
  }
  //polyData
  transformFilter->GetOutput()->GetCellData()->SetScalars(cellData);

  auto lookTable = vtkSmartPointer<vtkLookupTable>::New();
  lookTable->SetTableRange(bandMin, bandMax * koeff / 100);
  lookTable->SetHueRange(0, 0);
  lookTable->SetSaturationRange(0, 1);
  lookTable->SetValueRange(1, 1);
  lookTable->Build();

  auto scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
  scalarBar->SetLookupTable(lookTable);
  scalarBar->SetTitle("value");
  scalarBar->SetNumberOfLabels(Nband);

  /// Облучаемый объект - mapper
  auto bliskMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  bliskMapper->SetInputConnection(transformFilter->GetOutputPort());
  //bliskMapper->SetInputData(polyData);

  /// Облучаемый объект - actor
  auto bliskActor = vtkSmartPointer<vtkActor>::New();
  bliskActor->SetMapper(bliskMapper);
  bliskActor->GetProperty()->EdgeVisibilityOff();
  //09.12
  bliskActor->GetProperty()->SetAmbient(1);
  bliskActor->GetProperty()->SetDiffuse(0);
  bliskActor->GetProperty()->SetSpecular(0);

  auto endCamera = vtkSmartPointer<vtkCamera>::New();
  endCamera->SetViewUp(0, 1, 0);
  endCamera->SetPosition(-1, 0, 0.5);
  endCamera->SetFocalPoint(0, 0, 0);

  /// Visualize
  auto renderer = vtkSmartPointer<vtkRenderer>::New();
  auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(renderWindow);
  renderer->AddActor(bliskActor);
  renderer->AddActor2D(scalarBar);
  renderer->SetActiveCamera(endCamera);
  renderer->ResetCamera();
  //  renderer->SetBackground(1, 1, 1);
  renderWindow->SetSize(640, 480);
  renderer->ResetCameraClippingRange();
  renderWindow->Render();
  interactor->Initialize();

  ///Video mode
  if (VIDEOMODE)
  {
    auto callback = vtkSmartPointer<vtkTimerCallback2>::New();

    callback->actor = bliskActor;
    callback->mapper = bliskMapper;
    //callback->polyData = polyData;
    callback->polyData = transformFilter->GetOutput();
    callback->cellData = cellData;
    callback->colors = &interCells;

    callback->setAngleA(videoAngleA);
    callback->setAngleV(videoAngleV);
    callback->setTickLim(angleSaveLim);
    callback->setBandPiece(bandPiece);
    callback->setInitTick(kInit);
    callback->setNBand(Nband);

    interactor->AddObserver(vtkCommand::TimerEvent, callback);
    int timerId = interactor->CreateRepeatingTimer(videoTick);
    callback->timerId = timerId;

    system("pause");

    ///Video mode end
  }


  // Set the custom stype to use for interaction.
  vtkSmartPointer<MouseInteractorStyle> style =
    vtkSmartPointer<MouseInteractorStyle>::New();
  style->SetDefaultRenderer(renderer);
  style->Data = transformFilter->GetOutput();
  style->colors = &(*(interCells.end() - 1));

  interactor->SetInteractorStyle(style);

  /// Start the interaction and timer
  interactor->Start();
}

//___________________________________________________________
//____________________________main___________________________
int main(int argc, char* argv[])
{
  setlocale(LC_ALL, "Russian");
  //	SetConsoleOutputCP(1251);

  bool isCalculated = false;
  bool isWriting = false;
  std::string fileName;

  std::cout << "Считать записанные результаты (1) или посчитать (0)? ";
  std::cin >> isCalculated;

  vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
  reader->SetFileName("disc.stl");
  reader->Update();

  int numCells = reader->GetOutput()->GetNumberOfCells();


  int angleLim,
    angleSaveLim,
    angleTempLim;

  double angleAxis,    /// угол поворота относительно оси диска
    angleVert,    /// угол поворота относительно вертикали
    turnVert,    /// количество поворотов по вертикали (gamma)
    turnAxis,    /// количество поворотов по оси (alpha)
    koefAxisToVert;    /// коэффициент turnAxis / turnVert - 
       /// во сколько раз быстрее вокруг оси, чем отн-но вертикали
  int videoAngleV,
    videoAngleA;

  if (!isCalculated)
  {
    cout << "Введите число оборотов диска по вертикали: ";
    cin >> turnVert;
    cout << "Введите число оборотов диска по оси диска (положительное число)\n или коэффициент числа оборотов по оси диска к числу оборотов по вертикали \n (отрицательное число): ";
    cin >> turnAxis;
    cout << "Введите значение угла, через которое результаты будут сохраняться: ";
    cin >> videoAngleA;

    if (turnAxis < 0)
    {
      koefAxisToVert = -turnAxis;
      turnAxis = koefAxisToVert * turnVert;      // коэф-т, во сколько раз быстрее вокруг оси, чем отн-но вертикали
    }
    else
    {
      koefAxisToVert = turnAxis / turnVert; // коэф-т, во сколько раз быстрее вокруг оси, чем отн-но вертикали
    }

    if (turnAxis < 10)
      angleAxis = 0.5;
    else if (turnAxis < 100)
      angleAxis = 721. / 360;
    else
      angleAxis = 1081. / 360;
    angleVert = angleAxis / koefAxisToVert;

    angleLim = (int)(turnAxis * 360 / angleAxis);

    angleTempLim = (int)videoAngleA / angleAxis;
    videoAngleV = (int)videoAngleA / koefAxisToVert;
    angleSaveLim = angleLim / angleTempLim;
    if (angleLim % angleTempLim) angleSaveLim++;

  }
  else
  {
    std::cout << "Введите имя файла: ";
    std::cin >> fileName;

    std::ifstream resIn(fileName);
    if (!resIn)
    {
      std::cout << "Ошибка чтения файла" << std::endl;
      system("pause");
      return -1;
    }

    resIn >> angleSaveLim >> angleTempLim >> videoAngleA >> videoAngleV >> turnVert >> turnAxis >> koefAxisToVert;

    resIn.close();
  }

  std::vector <std::vector<double> > interCells(angleSaveLim, std::vector<double>(numCells));
  /// Данные об облученных ячейках: "доля" облучения (прореженная для сохранения)       

  double center[3], bounds[6], dimensions[3];

  auto initTransform = vtkSmartPointer<vtkTransform>::New();
  auto initTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  initTransformFilter->SetTransform(initTransform);
  initTransformFilter->SetInputConnection(reader->GetOutputPort());

  reader->GetOutput()->GetPoints()->GetBounds(bounds);
  dimensions[0] = std::abs(bounds[1] - bounds[0]);
  dimensions[1] = std::abs(bounds[3] - bounds[2]);
  dimensions[2] = std::abs(bounds[5] - bounds[4]);

  if (dimensions[0] < dimensions[1] && dimensions[0] < dimensions[2])
  {
    initTransform->RotateY(-90);
  }
  else if (dimensions[1] < dimensions[0] && dimensions[1] < dimensions[2])
  {
    initTransform->RotateX(90);
  }
  else
  {
    initTransform->Identity();
  }
  initTransformFilter->Update();

  auto initTransformTranslate = vtkSmartPointer<vtkTransform>::New();
  auto initTransformTranslateFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  initTransformTranslateFilter->SetTransform(initTransformTranslate);
  initTransformTranslateFilter->SetInputConnection(initTransformFilter->GetOutputPort());

  initTransformFilter->GetOutput()->GetCenter(center);
  initTransformTranslate->Translate(-center[0], -center[1], -center[2]);
  initTransformTranslateFilter->Update();

  initTransformTranslateFilter->GetOutput()->GetCenter(center);
  /*  std::cout << "Center of data is: "
              << center[0] << ", "
              << center[1] << ", "
              << center[2] << std::endl;*/
  initTransformTranslateFilter->GetOutput()->GetPoints()->GetBounds(bounds);
  /*  std::cout << "Bounds of data is: "
              << bounds[0] << ", "          //Xmin
              << bounds[1] << ", "          //Xmax
              << bounds[2] << ", "          //Ymin
              << bounds[3] << ", "          //Ymax
              << bounds[4] << ", "          //Zmin
              << bounds[5] << std::endl;    //Zmax*/

  auto normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
  normalGenerator->SetInputConnection(initTransformTranslateFilter->GetOutputPort());
  normalGenerator->ComputePointNormalsOff();
  normalGenerator->ComputeCellNormalsOn();
  normalGenerator->Update();

  auto regTransform = vtkSmartPointer<vtkTransform>::New();
  auto transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformFilter->SetTransform(regTransform);
  transformFilter->SetInputConnection(normalGenerator->GetOutputPort());
  transformFilter->Update();



  if (!isCalculated)
  {
    calculate(angleLim, angleTempLim, numCells, angleVert, angleAxis,
      interCells, transformFilter, regTransform);

    /// Сохранение результатов
    std::cout << "Сохранить результаты в файл?  ";
    std::cin >> isWriting;

    if (isWriting)
    {
      std::cout << "Введите имя файла: ";
      std::cin >> fileName;
      std::cout << "Запись в файл...  ";
      unsigned int start = GetTickCount();
      std::ofstream resOut(fileName);
      if (!resOut)
        std::cout << "Ошибка создания файла\n";
      else
      {
        resOut << angleSaveLim << ' ' << angleTempLim << ' ' << videoAngleA << ' ' << videoAngleV << ' '
          << turnVert << ' ' << turnAxis << ' ' << koefAxisToVert << std::endl;
        for (int k = 0; k < angleSaveLim; k++)
        {
          for (int i = 0; i < numCells; i++)
          {
            if (interCells[k][i] > 0)
              resOut << k << ' ' << i << ' ' << std::setprecision(2) << std::fixed << interCells[k][i] << std::endl;
          }
        }
        resOut.close();
      }
      unsigned int end = GetTickCount();
      std::cout << "Запись завершена" << std::endl;
      std::cout << "Длительность записи " << end - start << " ms" << std::endl;
    }
  }  // конец isCalculated
  else
  {
    readFile(fileName, numCells, interCells);
  }
  /// ========================= Отрисовка ==========================
  /// шкала - 8 градаций
  bool again = false;
  do
  {
    std::cout << "Подготовка к визуализации...";
    int Nband = 8;
    double bandMax = 0,//interCells[angleLim - 1][0], 
      bandMin = 0,//interCells[0][0], 
      bandPiece;

    for (auto& item : interCells[angleSaveLim - 1])
    {
      //  if (bandMin > item && item > 0) bandMin = item;
      if (bandMax < item)
      {
        bandMax = item;
      }
    }

    double koeff = 10;
    std::cout << "Введите количество градаций красного: ";
    std::cin >> Nband;
    std::cout << "Долю покрытия до какой величины (в % от максимума) разрешить шкалой? ";
    std::cin >> koeff;

    bandPiece = (bandMax - bandMin) / Nband * koeff / 100;

    visualize(bandMin, bandMax, bandPiece, Nband, koeff, numCells, interCells, transformFilter,
      angleSaveLim, videoAngleA, videoAngleV, angleTempLim);

    std::cout << "Хотите отрисовать результат еще раз? ";
    std::cin >> again;

  } while (again);

  if (!isCalculated && !isWriting)
  {
    std::cout << "Сохранить результаты в файл?  ";
    std::cin >> isWriting;

    if (isWriting)
    {
      std::cout << "Введите имя файла: ";
      std::cin >> fileName;
      std::cout << "Запись в файл...  ";
      unsigned int start = GetTickCount();
      std::ofstream resOut(fileName);
      if (!resOut)
        std::cout << "Ошибка создания файла\n";
      else
      {
        resOut << angleSaveLim << ' ' << angleTempLim << ' ' << videoAngleA << ' ' << videoAngleV << ' '
          << turnVert << ' ' << turnAxis << ' ' << koefAxisToVert << std::endl;
        for (int k = 0; k < angleSaveLim; k++)
        {
          for (int i = 0; i < numCells; i++)
          {
            if (interCells[k][i] > 0)
              resOut << k << ' ' << i << ' ' << std::setprecision(2) << std::fixed << interCells[k][i] << std::endl;
          }
        }
        resOut.close();
      }
      unsigned int end = GetTickCount();
      std::cout << "Запись завершена" << std::endl;
      std::cout << "Длительность записи " << end - start << " ms" << std::endl;
    }
  }

  return EXIT_SUCCESS;
}