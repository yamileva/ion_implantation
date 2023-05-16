cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

if(POLICY CMP0020)
  cmake_policy(SET CMP0020 NEW)
  cmake_policy(SET CMP0071 NEW)
endif()

project(TestModelFive)

find_package(VTK REQUIRED)

if (NOT VTK_FOUND)
  message("Skipping TestModelFive: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

if(vtkGUISupportQt_LOADED OR TARGET VTK::GUISupportQt)
  if (VTK_VERSION VERSION_LESS "8.90.0")
    # Instruct CMake to run moc automatically when needed.
    set(CMAKE_AUTOMOC ON)
  else()
     # Instruct CMake to run moc automatically when needed.
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
  endif()
  # We have ui files, this will bring in the macro: qt5_wrap_ui
  find_package(Qt5Widgets REQUIRED QUIET)
else()
  message(STATUS "GuiSupportQt not found.")
  return()
endif()

file(GLOB QT_WRAP *.h)
file(GLOB CXX_FILES *.cxx)

if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(TestModelFive MACOSX_BUNDLE ${CXX_FILES} ${QT_WRAP})
  qt5_use_modules(TestModelFive Core Gui)
  target_link_libraries(TestModelFive ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(TestModelFive MACOSX_BUNDLE ${CXX_FILES} ${QT_WRAP} )
  qt5_use_modules(TestModelFive Core Gui)
  target_link_libraries(TestModelFive ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS TestModelFive
    MODULES ${VTK_LIBRARIES}
    )
endif () 
