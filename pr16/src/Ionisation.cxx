#include <map>
#include <windows.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkLookupTable.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>  //09.12
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <vtksys/SystemTools.hxx>

#include "vtkTimerCallback2.h"
#include "Ionisation.h"

vtkStandardNewMacro(MouseInteractorStyle);


Ionisation::Ionisation() : QWidget()
{
  // значения по умолчанию
  turnVert = 80;
  koefAxisToVert = 1.2;
  turnAxis = turnVert * koefAxisToVert;

  angleAxis = 721. / 360.;
  videoAngleA = 50;

  vel = 1;
  velIon = 2 / 70 * 60;
  
  vizBandMax = 100;
  vizBandMin = 0;
  Nband = 8;
  
  //интерфейс
  resize(400, 400);
  buttonSTL = new QPushButton("Открыть STL файл");
  labelSTL = new QLabel("Файл не загружен");

  radioCalc = new QRadioButton("Выполнить новый расчет");
  radioRead = new QRadioButton("Отобразить готовые результаты");
  radioCalc->setChecked(true);
  radioCalc->setEnabled(false);
  radioRead->setEnabled(false);
  QGroupBox* groupIsCalc = new QGroupBox;
  QHBoxLayout* layIsCalc = new QHBoxLayout;

  layIsCalc->addWidget(radioCalc);
  layIsCalc->addWidget(radioRead);
  groupIsCalc->setLayout(layIsCalc);
  buttonRead = new QPushButton("Открыть файл результата");
  buttonRead->setEnabled(false);
  labelIsReaded = new QLabel("");
  progressReadBar = new QProgressBar();
  progressReadBar->setVisible(false);
  progressReadBar->setRange(0, 100);


  labelVel = new QLabel("Число оборотов диска по вертикали в минуту: ");
  labelStep = new QLabel("Шаг поворота (в градусах): ");
  labelRayCoeff = new QLabel("Коэффициент разрешения источника \n(1 - 10x40 точек): ");
  labelTurnVert = new QLabel("Число оборотов диска по вертикали: ");
  radioTurnAxis = new QRadioButton("Число оборотов диска по оси диска: ");
  radioTurnKoef = new QRadioButton("Отношение числа оборотов по оси диска к числу оборотов по вертикали: ");
  labelVideoAngleA = new QLabel("Период сохранения результатов (угол в градусах): ");
  
  radioTurnAxis->setChecked(true);
  QGroupBox* groupInput = new QGroupBox;
  QVBoxLayout* layInput = new QVBoxLayout;
  layInput->addWidget(radioTurnAxis);
  layInput->addWidget(radioTurnKoef);
  groupInput->setLayout(layInput);

  boxVel = new QDoubleSpinBox;
  boxVel->setValue(vel);
  boxVel->setRange(0.01, 10);
  boxStep = new QDoubleSpinBox;
  boxStep->setValue(angleAxis);
  boxStep->setRange(0.01, 10);
  boxRayCoeff = new QSpinBox;
  boxRayCoeff->setValue(RayCoeff);
  boxRayCoeff->setRange(1, 8);
  boxTurnVert = new QDoubleSpinBox;
  boxTurnVert->setValue(turnVert);
  boxTurnVert->setRange(1, 500);
  boxTurnAxis = new QDoubleSpinBox;
  boxTurnAxis->setValue(turnAxis);
  boxTurnAxis->setRange(1, 500);
  boxTurnKoef = new QDoubleSpinBox;
  boxTurnKoef->setValue(koefAxisToVert);
  boxTurnKoef->setRange(0.1, 50);
  boxTurnKoef->setEnabled(false);
  boxVideoAngleA = new QSpinBox;
  boxVideoAngleA->setValue(videoAngleA);
  boxVideoAngleA->setRange(10, 90);
  

  buttonCalc = new QPushButton("Рассчитать");
  buttonCalc->setEnabled(false);
  progressBar = new QProgressBar();
  progressBar->setVisible(false);
  progressBar->setRange(0, 100);
  labelCalcStatus = new QLabel("");

  buttonWrite = new QPushButton("Сохранить результат в файл");
  buttonWrite->setEnabled(false);
  labelIsWrited = new QLabel("");
  progressWriteBar = new QProgressBar();
  progressWriteBar->setVisible(false);
  progressWriteBar->setRange(0, 100);


  labelNband = new QLabel("Количество градаций шкалы: ");
  labelBandMax = new QLabel("Верхняя граница шкалы (% от максимума): ");
  labelBandMin = new QLabel("Нижняя граница шкалы (% от максимума): ");
  
  boxNband = new QSpinBox;
  boxNband->setValue(Nband);
  boxNband->setRange(4, 20);
  boxBandMax = new QDoubleSpinBox;
  boxBandMax->setValue(vizBandMax);
  boxBandMax->setRange(0, 100);
  boxBandMin = new QDoubleSpinBox;
  boxBandMin->setValue(vizBandMin);
  boxBandMin->setRange(0, 100);
  
  buttonVis = new QPushButton("Отобразить результаты");
  buttonVis->setEnabled(false);
  checkVideo = new QCheckBox("Анимация");
  checkVideo->setChecked(false);
  //labelAnim = new QLabel("Анимация");

  QGridLayout* allLayout = new QGridLayout;
  int i = 0;
  allLayout->addWidget(buttonSTL, i, 0);
  allLayout->addWidget(labelSTL, i, 1);
  i++;
  allLayout->addWidget(groupIsCalc, i, 0, 1, 2);
  i++;
  allLayout->addWidget(buttonRead, i, 0);
  allLayout->addWidget(labelIsReaded, i, 1);
  allLayout->addWidget(progressReadBar, i, 1);
  i++;
  allLayout->addWidget(labelVel, i, 0, 1, 2);
  allLayout->addWidget(boxVel, i, 2);
  i++;
  allLayout->addWidget(labelStep, i, 0, 1, 2);
  allLayout->addWidget(boxStep, i, 1);
  i++;
  allLayout->addWidget(labelRayCoeff, i, 0, 1, 2);
  allLayout->addWidget(boxRayCoeff, i, 1);
  i++;
  allLayout->addWidget(labelTurnVert, i, 0, 1, 2);
  allLayout->addWidget(boxTurnVert, i, 2);
  i++;
  allLayout->addWidget(groupInput, i, 0, 2, 2);
  allLayout->addWidget(boxTurnAxis, i, 2);
  i++; 
  allLayout->addWidget(boxTurnKoef, i, 2);
  i++;
  allLayout->addWidget(labelVideoAngleA, i, 0, 1, 2);
  allLayout->addWidget(boxVideoAngleA, i, 2);
  i++;
  allLayout->addWidget(buttonCalc, i, 0);
  allLayout->addWidget(progressBar, i, 1, 1, 2);
  i++;
  allLayout->addWidget(labelCalcStatus, i, 1);
  i++;
  allLayout->addWidget(buttonWrite, i, 0);
  allLayout->addWidget(labelIsWrited, i, 1);
  allLayout->addWidget(progressWriteBar, i, 1);
  i++;
  allLayout->addWidget(labelNband, i, 0, 1, 2);
  allLayout->addWidget(boxNband, i, 2);
  i++;
  allLayout->addWidget(labelBandMax, i, 0, 1, 2);
  allLayout->addWidget(boxBandMax, i, 2);
  i++;
  allLayout->addWidget(labelBandMin, i, 0, 1, 2);
  allLayout->addWidget(boxBandMin, i, 2);
  i++;
  allLayout->addWidget(buttonVis, i, 0);
  //allLayout->addWidget(labelAnim, i, 1);
  allLayout->addWidget(checkVideo, i, 1);
  
  setLayout(allLayout);

  
  // connections

  QObject::connect(buttonSTL, SIGNAL(clicked()), SLOT(slotButtonSTLClicked()));
  QObject::connect(buttonRead, SIGNAL(clicked()), SLOT(slotButtonReadClicked()));
  QObject::connect(buttonCalc, SIGNAL(clicked()), SLOT(slotButtonCalcClicked()));
  QObject::connect(buttonWrite, SIGNAL(clicked()), SLOT(slotButtonWriteClicked()));
  QObject::connect(buttonVis, SIGNAL(clicked()), SLOT(slotButtonVisClicked()));

  QObject::connect(radioCalc, SIGNAL(clicked()), SLOT(slotRadioCalcChecked()));
  QObject::connect(radioRead, SIGNAL(clicked()), SLOT(slotRadioReadChecked()));
  QObject::connect(radioTurnAxis, SIGNAL(clicked()), SLOT(slotRadioTurnAxisChecked()));
  QObject::connect(radioTurnKoef, SIGNAL(clicked()), SLOT(slotRadioTurnKoefChecked()));
  QObject::connect(checkVideo, SIGNAL(clicked()), SLOT(slotCheckVideoChecked()));
 

  regTransform = vtkSmartPointer<vtkTransform>::New();
  transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformFilter->SetTransform(regTransform);

}

void Ionisation::initGeometry()
{
  QByteArray ba = STLName.toLocal8Bit();
  const char* c_str2 = ba.data();
  vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
  reader->SetFileName(c_str2);
  reader->Update();

  numCells = reader->GetOutput()->GetNumberOfCells();

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

  transformFilter->SetInputConnection(normalGenerator->GetOutputPort());
  transformFilter->Update();
  isSTLOpened = true;


}

void Ionisation::visualize()
{
  auto cellData = vtkSmartPointer<vtkUnsignedCharArray>::New();
  cellData->SetNumberOfComponents(3);
  cellData->SetNumberOfTuples(numCells);

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


  bandPiece = bandMax * (vizBandMax - vizBandMin) / Nband / 100;

  bandMin = bandMax * vizBandMin / 100;
  
  int kInit;
  if (isAnimating)
    kInit = 0;
  else
    kInit = angleSaveLim - 1;

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
  transformFilter->GetOutput()->GetCellData()->SetScalars(cellData);

  auto lookTable = vtkSmartPointer<vtkLookupTable>::New();
  lookTable->SetTableRange(bandMin, bandMax * vizBandMax / 100);
  lookTable->SetHueRange(0, 0);
  lookTable->SetSaturationRange(0, 1);
  lookTable->SetValueRange(1, 1);
  lookTable->Build();

  auto scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
  scalarBar->SetLookupTable(lookTable);
  scalarBar->SetTitle("x10^17");
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

  if (isAnimating)
  {
    auto callback = vtkSmartPointer<vtkTimerCallback2>::New();

    callback->actor = bliskActor;
    callback->mapper = bliskMapper;
    //callback->polyData = polyData;
    callback->polyData = transformFilter->GetOutput();
    callback->cellData = cellData;
    callback->colors = &interCells;

    callback->setAngleA(videoAngleA);
    callback->setTickLim(angleSaveLim);
    callback->setBandPiece(bandPiece);
    callback->setInitTick(kInit);
    callback->setNBand(Nband);

    interactor->AddObserver(vtkCommand::TimerEvent, callback);
    int timerId = interactor->CreateRepeatingTimer(angleTempLim);
    callback->timerId = timerId;
  }
  ///Video mode end
  vtkSmartPointer<MouseInteractorStyle> style =
    vtkSmartPointer<MouseInteractorStyle>::New();
  style->SetDefaultRenderer(renderer);
  style->SetCurrentRenderer(renderer);
  style->Data = transformFilter->GetOutput();
  style->colors = &(*(interCells.end() - 1));

  interactor->SetInteractorStyle(style);

  /// Start the interaction and timer

  interactor->Start();

}

void Ionisation::calculate()
{
  int rayCountX = 11 * RayCoeff;
  int rayCountY = 41 * RayCoeff;
  int rayCount = rayCountX * rayCountY * 2;

  double rayWidthTotal = 100, rayHeightTotal = 400;
  double rayWidth = rayWidthTotal / (rayCountX - 1),
    rayHeight = rayHeightTotal / (rayCountY - 1);

  double rayOrigin[3] = { -rayWidthTotal / 2, 150, 600 };
  double axis[3] = { 0.0, 0.0, 1.0 },
    temp[3];

  std::vector <std::vector<double> >  rayStart(rayCount, std::vector<double>(3)), 
                                      rayEnd(rayCount, std::vector<double>(3));

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
  
  double xyz[3], t, pcoords[3], attack, cellNormal[3], singleRayStart[3], singleRayEnd[3];
  double angleSaveMoment = (angleLim - 1) % angleTempLim;
  int subId, hit = 0;
  vtkIdType cellId = -1;
  int k, kTemp, kSave;
  double deltaT = angleVert / vel / 360;  // в минутах
  progressBar->setValue(0);


  for (k = 0, kTemp = 1, kSave = 0; k < angleLim; k++, kTemp++)  /// цикл по "времени"
  {
    cellLocator->BuildLocator();
    hit = 0;

    for (int i = 0; i < rayCount; i++)  /// цикл по точкам источника
    {
      singleRayStart[0] = rayStart[i][0];
      singleRayStart[1] = rayStart[i][1];
      singleRayStart[2] = rayStart[i][2];
      singleRayEnd[0] = rayEnd[i][0];
      singleRayEnd[1] = rayEnd[i][1];
      singleRayEnd[2] = rayEnd[i][2];
      hit = cellLocator->IntersectWithLine(singleRayStart, singleRayEnd, 0.0001, t, xyz, pcoords, subId, cellId);

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
        interCellsTemp.at(kTemp)[item.first] += item.second.second / item.second.first * deltaT * velIon;
    }
    else
    {
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

    if ((k + 1) % (angleLim / 100) == 0)
    {
      progressBar->setValue((k + 1) * 100 / angleLim);
    }

    /// поворот диска вокруг своей оси
    regTransform->RotateWXYZ(angleAxis, axis);
    transformFilter->Update();



  } // конец цикла по k


}

int Ionisation::readFile()
{
  QByteArray ba = fileName.toLocal8Bit();
  const char* c_str2 = ba.data();
  
  std::ifstream resIn(c_str2);
  if (!resIn)
  {
    return -1;
  }
  resIn >> angleSaveLim >> angleTempLim >> videoAngleA 
        >> turnVert >> turnAxis >> koefAxisToVert;

  interCells = std::vector <std::vector<double> >(angleSaveLim, std::vector<double>(numCells, 0));

  int k, i;
  double num;
  progressReadBar->setValue(0);
  while (!resIn.eof())
  {
    resIn >> k >> i >> num;
    interCells[k][i] = num;
    if ((k + 1) % (angleSaveLim / 100) == 0)
    {
      progressReadBar->setValue((k + 1) * 100 / angleSaveLim);
    }
  }
  resIn.close();

  return 0;
}

int Ionisation::writeFile()
{
  QByteArray ba = fileNameW.toLocal8Bit();
  const char* c_str2 = ba.data();
  
  std::ofstream resOut(c_str2);
  if (!resOut)
    return -1;
  else
  {
    progressWriteBar->setValue(0);
    resOut << angleSaveLim << ' ' << angleTempLim << ' ' << videoAngleA << ' '
      << turnVert << ' ' << turnAxis << ' ' << koefAxisToVert << std::endl;
    for (int k = 0; k < angleSaveLim; k++)
    {
      for (int i = 0; i < numCells; i++)
      {
        if (interCells[k][i] > 0)
          resOut << k << ' ' << i << ' ' << std::setprecision(2) << std::fixed << interCells[k][i] << std::endl;
      }
      if ((k + 1) % (angleSaveLim / 100) == 0)
      {
        progressWriteBar->setValue((k + 1) * 100 / angleSaveLim);
      }
    }
    resOut.close();
  }
  return 0;
}

void Ionisation::slotButtonSTLClicked()
{
  STLName = QFileDialog::getOpenFileName(0, "Загрузить файл STL", "", "*.stl");
  if (STLName.isEmpty())
  {
    labelSTL->setText("Ошибка выбора файла");
    buttonVis->setEnabled(false);
    return;
  }
  
  labelSTL->setText(STLName);
  isSTLOpened = true;
  initGeometry();
  radioCalc->setEnabled(true);
  radioRead->setEnabled(true);
  buttonCalc->setEnabled(true);

}

void Ionisation::slotButtonReadClicked()
{
  fileName = QFileDialog::getOpenFileName(0, "Загрузить файл результатов", "", "*.*");
  if (fileName.isEmpty())
  {
    labelIsReaded->setText("Ошибка выбора файла");
    buttonVis->setEnabled(false);
    return;
  }

  labelIsReaded->setText("");
  labelIsReaded->setVisible(false);
  progressReadBar->reset();
  progressReadBar->setVisible(true);

  int res = readFile();
  progressReadBar->setVisible(false);
  labelIsReaded->setVisible(true);
  if (res == -1)
  {
    labelIsReaded->setText("Ошибка чтения файла");
    buttonVis->setEnabled(false);
  }
  else
  {
    buttonVis->setEnabled(true);
    labelIsReaded->setText(fileName);
    boxTurnVert->setValue(turnVert);
    boxTurnAxis->setValue(turnAxis);
    boxTurnKoef->setValue(koefAxisToVert);
    boxVideoAngleA->setValue(videoAngleA);
    
  }
  

}

void Ionisation::slotButtonCalcClicked()
{
  labelCalcStatus->setText("Входные данные изменены");
  turnVert = boxTurnVert->value();
  vel = boxVel->value();

  if (radioTurnAxis->isChecked())
  {
    turnAxis = boxTurnAxis->value();
    koefAxisToVert = turnAxis / turnVert;
  }
  else
  {
    koefAxisToVert = boxTurnKoef->value();
    turnAxis = koefAxisToVert * turnVert;
  }
  videoAngleA = boxVideoAngleA->value();


  angleAxis = boxStep->value();
  angleVert = angleAxis / koefAxisToVert;
  RayCoeff = boxRayCoeff->value();

  angleLim = (int)(turnAxis * 360 / angleAxis);
  
  angleTempLim = (int)videoAngleA / angleAxis;
  angleSaveLim = angleLim / angleTempLim;
  if (angleLim % angleTempLim) angleSaveLim++;
  labelCalcStatus->setText("");
  std::cout << "AngleLim " << angleLim << endl;
  std::cout << "angleTempLim " << angleTempLim << endl;
  std::cout << "angleSaveLim " << angleSaveLim << endl;

  if (angleSaveLim > 0)
  {
    if (interCells.size() != 0)
    {
      interCells.clear();
    }
    interCells = std::vector <std::vector<double> >(angleSaveLim, std::vector<double>(numCells, 0));
  }
  else
    return;
  
  progressBar->reset();
  progressBar->setVisible(true);
  calculate();
  progressBar->setVisible(false);
  labelCalcStatus->setText("Расчет окончен");
  buttonVis->setEnabled(true);
  buttonWrite->setEnabled(true);

}

void Ionisation::slotButtonWriteClicked()
{
  fileNameW = QFileDialog::getSaveFileName(0, "Сохранить файл результатов", "*.txt", "*.txt");
  if (fileNameW.isEmpty())
  {
    labelIsWrited->setText("Ошибка выбора файла");
    return;
  }

  labelIsWrited->setText("");
  labelIsWrited->setVisible(false);
  progressWriteBar->reset();
  progressWriteBar->setVisible(true);

  int res = writeFile();
  progressWriteBar->setVisible(false);
  labelIsWrited->setVisible(true);
  
  if (res == -1)
  {
    labelIsWrited->setText("Ошибка записи файла");
  }
  else
    labelIsWrited->setText(fileNameW);

}

void Ionisation::slotButtonVisClicked()
{
  Nband = boxNband->value();
  vizBandMax = boxBandMax->value();
  vizBandMin = boxBandMin->value();
  visualize();
}

void Ionisation::slotRadioCalcChecked()
{
  buttonRead->setEnabled(false);

}

void Ionisation::slotRadioReadChecked()
{
  buttonRead->setEnabled(true);
}

void Ionisation::slotRadioTurnAxisChecked()
{
  boxTurnAxis->setEnabled(true);
  boxTurnKoef->setEnabled(false);
}

void Ionisation::slotRadioTurnKoefChecked()
{
  boxTurnAxis->setEnabled(false);
  boxTurnKoef->setEnabled(true);
}

void Ionisation::slotCheckVideoChecked()
{
  isAnimating = checkVideo->isChecked();
}