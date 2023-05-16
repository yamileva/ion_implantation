#pragma once

#include <QtWidgets>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>


class Ionisation :
  public QWidget
{
  Q_OBJECT
  public:
    Ionisation();
    ~Ionisation() = default;
    void initGeometry();
    void visualize();
    int readFile();
    int writeFile();
    void calculate();


  public slots:
    void slotButtonSTLClicked();
    void slotButtonReadClicked();
    void slotButtonCalcClicked();
    void slotButtonWriteClicked();
    void slotButtonVisClicked();

    void slotRadioCalcChecked();
    void slotRadioReadChecked();

    void slotRadioTurnAxisChecked();
    void slotRadioTurnKoefChecked();

    void slotCheckVideoChecked();

  private:

    QPushButton* buttonSTL;
    QLabel* labelSTL;
    
    QRadioButton* radioCalc;
    QRadioButton* radioRead;
    QPushButton* buttonRead;
    QLabel* labelIsReaded;
    QProgressBar* progressReadBar;

    QLabel* labelVel;
    QLabel* labelStep;
    QLabel* labelRayCoeff;

    QLabel* labelTurnVert;
    QRadioButton* radioTurnAxis;
    QRadioButton* radioTurnKoef;
    QLabel* labelVideoAngleA;

    QDoubleSpinBox* boxVel; 
    QDoubleSpinBox* boxStep;
    QSpinBox* boxRayCoeff;
    QDoubleSpinBox* boxTurnVert;
    QDoubleSpinBox* boxTurnAxis;
    QDoubleSpinBox* boxTurnKoef;
    QSpinBox* boxVideoAngleA;

    QPushButton* buttonCalc;
    QProgressBar* progressBar;
    QLabel* labelCalcStatus;

    QPushButton* buttonWrite;
    QLabel* labelIsWrited;
    QProgressBar* progressWriteBar;


    QLabel* labelNband;
    QLabel* labelBandMax;
    QLabel* labelBandMin;

    QSpinBox* boxNband;
    QDoubleSpinBox* boxBandMax;
    QDoubleSpinBox* boxBandMin;

    QCheckBox* checkVideo;
    QPushButton* buttonVis;
    QLabel* labelAnim;


    bool isSTLOpened = false;
    bool isCalculated = false;
    bool isWriting = false;
    bool isAnimating = false;

    QString STLName;
    QString fileName;
	  QString fileNameW;

    int angleLim = 0,
      angleSaveLim = 0,
      angleTempLim = 0,
      numCells = 0,
      RayCoeff = 2;

    double angleAxis,    /// угол поворота относительно оси диска
      angleVert,    /// угол поворота относительно вертикали
      turnVert,    /// количество поворотов по вертикали (gamma)
      turnAxis,    /// количество поворотов по оси (alpha)
      koefAxisToVert,    /// коэффициент turnAxis / turnVert - 
         /// во сколько раз быстрее вокруг оси, чем отн-но вертикали
      vel,  /// количество оборотов вокруг вертикальной оси в минуту
      velIon;  /// доза облучения на 1 см^2 в минуту
    double videoAngleA;

    double vizBandMax,
        vizBandMin;

    int Nband;
    
    std::vector <std::vector<double> > interCells;

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter;
    vtkSmartPointer<vtkTransform> regTransform;


};

