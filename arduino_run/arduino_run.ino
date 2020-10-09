#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <MemoryFree.h>
#include <pgmStrToRAM.h>
#define grid 13
#define gridheight 13



void setup() {
  Serial.begin(115200);
}
void loop() {
  Serial.println(F("Starting..."));
  
  double u[grid][gridheight+1], un[grid][gridheight+1], uc[grid][gridheight+1], un2[grid][gridheight+1];
  double v[grid+1][gridheight], vn[grid+1][gridheight], vc[grid+1][gridheight], vn2[grid+1][gridheight];
  double p[grid+1][gridheight+1], pn[grid+1][gridheight+1], pc[grid+1][gridheight+1], pn2[grid+1][gridheight+1];
  double m[grid+1][gridheight+1];
  int i, j, stepVal;
  double dx, dy, dt, tau, delta, errorVal, Re;
  int maxstep = 10000;
  double residual[maxstep], umomresidual[maxstep], vmomresidual[maxstep];
  stepVal = 1;
  dx = 1.0/(grid-1);
  dy = dx;
  //dt = 0.0005;
  dt=0.0005;
  delta = 1.0;
  errorVal = 1.0;
  Re = 3200.0;
  unsigned int start, endVal;
  double cpu_time_used;
  double CFL1, CFL2;
  CFL1 = dt/dx;

  int precision = 6;

  char CFL1buf[20];
  char CFL2buf[20];
  dtostrf(CFL1, precision+3, precision, CFL1buf);
  dtostrf(CFL2, precision+3, precision, CFL2buf);
  Serial.print(F("CFL numbers are: "));
  Serial.print(CFL1buf);
  Serial.print(F(", "));
  Serial.print(CFL2buf);
  Serial.println();

  double WU;
  // Initializing u
    for (i=0; i<=(grid-1); i++)
    {
      for (j=0; j<=(gridheight); j++)
      {
        u[i][j] = 0.0;
        u[i][gridheight] = 1.0;
        u[i][gridheight-1] = 1.0;
      }
    }
    
  // Initializing v
    for (i=0; i<=(grid); i++)
    {
      for (j=0; j<=(gridheight-1); j++)
      {
        v[i][j] = 0.0;
      }
    }
    
  // Initializing p
    for (i=0; i<=(grid); i++)
    {
      for (j=0; j<=(gridheight); j++)
      {
        p[i][j] = 1.0;
      }
    }
  
  start = millis();
  while (errorVal > 0.000001 && stepVal<=maxstep)
  {
    // Solve u-momentum
    for (i=1; i<=(grid-2); i++)
    {
      for (j=1; j<=(gridheight-1); j++)
      {
        un[i][j] = u[i][j] - dt*(  (u[i+1][j]*u[i+1][j]-u[i-1][j]*u[i-1][j])/2.0/dx 
              +0.25*( (u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i][j]+u[i][j-1])*(v[i+1][j-1]+v[i][j-1]) )/dy  )
                - dt/dx*(p[i+1][j]-p[i][j]) 
                  + dt*1.0/Re*( (u[i+1][j]-2.0*u[i][j]+u[i-1][j])/dx/dx +(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dy/dy );
      }
    }
    
    for (j=1; j<=(gridheight-1); j++)
    {
      un[0][j] = 0.0;
      un[grid-1][j] = 0.0;
    }
    
    for (i=0; i<=(grid-1); i++)
    {
      un[i][0] = -un[i][1];
      un[i][gridheight] = 2.0 - un[i][gridheight-1];
    }
    
    
    // Solves v-momentum
    for (i=1; i<=(grid-1); i++)
    {
      for (j=1; j<=(gridheight-2); j++)
      {
        vn[i][j] = v[i][j] - dt* ( 0.25*( (u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i-1][j]+u[i-1][j+1])*(v[i][j]+v[i-1][j]) )/dx 
              +(v[i][j+1]*v[i][j+1]-v[i][j-1]*v[i][j-1])/2.0/dy ) 
                - dt/dy*(p[i][j+1]-p[i][j]) 
                  + dt*1.0/Re*( (v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy/dy );
      }
    }
    
    for (j=1; j<=(gridheight-2); j++)
    {
      vn[0][j] = -vn[1][j];
      vn[grid][j] = -vn[grid-1][j];
    }   

    for (i=0; i<=(grid); i++)
    {
      vn[i][0] = 0.0;
      vn[i][gridheight-1] = 0.0;
    }   
  
    // Solves continuity equation
    for (i=1; i<=(grid-1); i++)
    {
      for (j=1; j<=(gridheight-1); j++)
      {
        pn[i][j] = p[i][j]-dt*delta*(  ( un[i][j]-un[i-1][j] )/dx + ( vn[i][j]-vn[i][j-1] ) /dy  );
      }
    }
    
    for (i=1; i<=(grid-1); i++)
    {
      pn[i][0] = pn[i][1];
      pn[i][gridheight] = pn[i][gridheight-1];
    }
    
    for (j=0; j<=(gridheight); j++)
    {
      pn[0][j] = pn[1][j];
      pn[grid][j] = pn[grid-1][j];
    }   
    
    // Displaying error
    errorVal = 0.0;
    
    for (i=1; i<=(grid-1); i++)
    {
      for (j=1; j<=(gridheight-1); j++)
      {
        m[i][j] = (  ( un[i][j]-un[i-1][j] )/dx + ( vn[i][j]-vn[i][j-1] )/dy  );
        errorVal = errorVal + fabs(m[i][j]);
      }
    }
    // residual[step] = log10(errorVal);
    
    if (stepVal%10 == 0 && (stepVal!=0))
    {
      char errorBuf[20];
      char stepBuf[20];
      int precision = 6;
      dtostrf(errorVal, precision+3, precision, errorBuf);
      dtostrf(stepVal, precision+3, precision, stepBuf);
      Serial.print(F("Error is: "));
      Serial.print(errorBuf);
      Serial.print(F(" for the step: "));
      Serial.print(stepBuf);
      Serial.println();
    }
    
    // Iterating u
    for (i=0; i<=(grid-1); i++)
    {
      for (j=0; j<=(gridheight); j++)
      {
        u[i][j] = un[i][j];
      }
    }
    
    // Iterating v
    for (i=0; i<=(grid); i++)
    {
      for (j=0; j<=(gridheight-1); j++)
      {
        v[i][j] = vn[i][j];
      }
    }
    
    // Iterating p
    for (i=0; i<=(grid); i++)
    {
      for (j=0; j<=(gridheight); j++)
      {
        p[i][j] = pn[i][j];
      }
    }

    stepVal = stepVal + 1;
    WU = WU + 1;
  }
  endVal = millis();
  cpu_time_used = ((double) (endVal - start)) / 16000000; //16MHz for MEGA2560
  
  for (i=0; i<=(grid-1); i++)
    {
      for (j=0; j<=(gridheight-1); j++)
      { 
          uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
                vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
                pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
      }
    }
    
  char cpuBuf[20];
  char wuBuf[20];
  dtostrf(cpu_time_used, precision+3, precision, cpuBuf);
  dtostrf(WU, precision+3, precision, wuBuf);
  Serial.print(F("CPU Time = "));
  Serial.print(cpuBuf);
  Serial.println();
  Serial.print(F("Workunits = "));
  Serial.print(wuBuf);
  Serial.println();
  
  
  // OUTPUT DATA
  Serial.println(F("VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\""));
  Serial.println(F("ZONE  F=POINT\n"));
  Serial.print(F("I="));
  Serial.print(grid, DEC);
  Serial.print(F(", J="));
  Serial.println(gridheight, DEC);
  

  Serial.println(F("XPos\t\tYPos\t\tUC\t\tVC\t\tPC"));
  for ( j = 0 ; j < (gridheight) ; j++ )
  {
    for ( i = 0 ; i < (grid) ; i++ )
    {
    double xpos, ypos;
    xpos = i*dx;
    ypos = j*dy;
    
    Serial.print(xpos, DEC);
    Serial.print("\t");
    Serial.print(ypos, DEC);
    Serial.print("\t");
    Serial.print(uc[i][j], DEC);
    Serial.print("\t");
    Serial.print(vc[i][j], DEC);
    Serial.print("\t");
    Serial.println(pc[i][j], DEC);
    }
  }
  
  // CENTRAL --U
    Serial.println(F("VARIABLES=\"U\",\"Y\""));
    Serial.println(F("ZONE F=POINT"));
    Serial.print(F("I="));
    Serial.println(gridheight, DEC);

    for ( j = 0 ; j < gridheight ; j++ )
    {
      double ypos;
      ypos = (double) j*dy;
  
      Serial.print((uc[(grid-1)/2][j]), DEC);
      Serial.print("\t");
      Serial.println(ypos, DEC);
    }
   
  // CENTRAL --V
    Serial.println(F("VARIABLES=\"V\",\"X\""));
    Serial.println(F("ZONE F=POINT"));
    Serial.print(F("I="));
    Serial.println(grid, DEC);

    for ( i = 0 ; i < grid ; i++ )
    {
      double xpos;
      xpos = (double) i*dx;
  
      Serial.print((vc[i][(gridheight-1)/2]), DEC);
      Serial.print("\t");
      Serial.println(xpos, DEC);
    }  
  // Residual
    // fprintf(fout4, "VARIABLES=\"Workunits\",\"Error\"\n");
    // fprintf(fout4, "ZONE F=POINT\n");
    // fprintf(fout4, "I=%d\n", maxstep );

    // for ( j = 0 ; j < maxstep ; j++ )
    // {
    // fprintf( fout4, "%d\t%5.18lf\n", j, residual[j] );
    // }
    
  // // U Momentum Residual
    // fprintf(fout5, "VARIABLES=\"Step\",\"U-Res\"\n");
    // fprintf(fout5, "ZONE F=POINT\n");
    // fprintf(fout5, "I=%d\n", maxstep );

    // for ( j = 0 ; j < maxstep ; j++ )
    // {
    // fprintf( fout5, "%d\t%5.18lf\n", j, umomresidual[j] );
    // }
    
  // // V Momentum Residual
    // fprintf(fout6, "VARIABLES=\"Step\",\"V-Res\"\n");
    // fprintf(fout6, "ZONE F=POINT\n");
    // fprintf(fout6, "I=%d\n", maxstep );

    // for ( j = 0 ; j < maxstep ; j++ )
    // {
    // fprintf( fout6, "%d\t%5.18lf\n", j, vmomresidual[j] );
    // }  
    
    Serial.println(F("Finished."));
    while(1);
}
