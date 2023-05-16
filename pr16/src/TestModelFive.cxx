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
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <vtksys/SystemTools.hxx>

#include "Ionisation.h"


int main(int argc, char* argv[])
{
  //setlocale(LC_ALL, "Russian");
  //	SetConsoleOutputCP(1251);
  //QTextCodec::setCodecForCStrings(QTextCodec::codecForName("UTF-8")); //изменения
  QTextCodec::setCodecForLocale(QTextCodec::codecForName("UTF-8")); //изменения
  //QTextCodec::setCodecForTr(QTextCodec::codecForName("UTF-8")); //изменения

  QApplication app(argc, argv);
  Ionisation ion;

  ion.show();

  return app.exec();
}