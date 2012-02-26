/*
 * =====================================================================================
 *
 *       Filename:  Visualization.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/13/2010 10:04:33 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#include "Global.h"

#ifndef VISUALIZATION_VTK
#define VISUALIZATION_VTK

#include "Setup.h"


/* 
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h" 
#include "vtkProperty.h" 
#include "vtkActor.h" 
#include "vtkRenderWindow.h" 
#include "vtkRenderer.h" 
#include "vtkRenderWindowInteractor.h"
#include "vtkWindowToImageFilter.h"
#include "vtkJPEGWriter.h"
*/


class Visualization {

  public:
    Visualization(Setup *setup) {
    
    
    
    };

    ~Visualization() {



    };

    int visualize() {

  // create sphere geometry 
  vtkSphereSource *sphere = vtkSphereSource::New(); 
  sphere->SetRadius(1.0); 
  sphere->SetThetaResolution(18); 
  sphere->SetPhiResolution(18);

  // map to graphics library 
  vtkPolyDataMapper *map = vtkPolyDataMapper::New(); 
  map->SetInput(sphere->GetOutput());

  // actor coordinates geometry, properties, transformation 
  vtkActor *aSphere = vtkActor::New(); 
  aSphere->SetMapper(map); 
  aSphere->GetProperty()->SetColor(0,0,1); // sphere color blue

  // a renderer and render window 
  vtkRenderer *ren1 = vtkRenderer::New(); 
  vtkRenderWindow *renWin = vtkRenderWindow::New(); 
  renWin->AddRenderer(ren1);
  renWin->SetOffScreenRendering(true); 

  // an interactor 
//  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New(); 
//  iren->SetRenderWindow(renWin);

  // add the actor to the scene 
  ren1->AddActor(aSphere); 
  ren1->SetBackground(1,1,1); // Background color white

  // render an image (lights and cameras are created automatically) 
  renWin->Render();


    vtkWindowToImageFilter * w2i  = vtkWindowToImageFilter::New(); 
    w2i->SetInput( renWin ); 
    vtkJPEGWriter * image  = vtkJPEGWriter::New(); 
    w2i->Modified(); 
    image->SetInput( w2i->GetOutput()); 
    image->SetFileName( "aaa.jpg" ); 
    image->SetQuality ( 100 ); 
    image->Write (); 


  // release memory and return 
  sphere->Delete(); 
  map->Delete(); 
  aSphere->Delete(); 
  ren1->Delete(); 
  renWin->Delete(); 
  //iren->Delete(); 
  return EXIT_SUCCESS; 

    }


};


#endif // VISUALIZATION_VTK
 * */
