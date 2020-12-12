/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _JELLO_H_
#define _JELLO_H_

#include "openGL-headers.h"
#include "pic.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define pi 3.141592653589793238462643383279 

// camera angles
extern double Theta;
extern double Phi;
extern double R;

// number of images saved to disk so far
extern int sprite;

// mouse control
extern int g_vMousePos[2];
extern int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

struct point 
{
   double x;
   double y;
   double z;
};

// these variables control what is displayed on the screen
extern int shear, bend, structural, pause, viewingMode, saveScreenToFile;

struct world
{
  char integrator[10]; // "RK4" or "Euler"
  double dt; // timestep, e.g.. 0.001
  int n; // display only every nth timepoint
  double kElastic; // Hook's elasticity coefficient for all springs except collision springs
  double dElastic; // Damping coefficient for all springs except collision springs
  double kCollision; // Hook's elasticity coefficient for collision springs
  double dCollision; // Damping coefficient collision springs
  double mass; // mass of each of the 512 control points, mass assumed to be equal for every control point
  int incPlanePresent; // Is the inclined plane present? 1 = YES, 0 = NO (always NO in this assignment)
  double a,b,c,d; // inclined plane has equation a * x + b * y + c * z + d = 0; if no inclined plane, these four fields are not used
  int resolution; // resolution for the 3d grid specifying the external force field; value of 0 means that there is no force field
  struct point * forceField; // pointer to the array of values of the force field
  struct point p[8][8][8]; // position of the 512 control points
  struct point v[8][8][8]; // velocities of the 512 control points
  char * fileName;
};

extern struct world jello;

// computes crossproduct of three vectors, which are given as points
// struct point vector1, vector2, dest
// result goes into dest
#define CROSSPRODUCTp(vector1,vector2,dest)\
  CROSSPRODUCT( (vector1).x, (vector1).y, (vector1).z,\
                (vector2).x, (vector2).y, (vector2).z,\
                (dest).x, (dest).y, (dest).z )

// computes crossproduct of three vectors, which are specified by floating-point coordinates
// double x1,y1,z1,x2,y2,z2,x,y,z
// result goes into x,y,z
#define CROSSPRODUCT(x1,y1,z1,x2,y2,z2,x,y,z)\
\
  x = (y1) * (z2) - (y2) * (z1);\
  y = (x2) * (z1) - (x1) * (z2);\
  z = (x1) * (y2) - (x2) * (y1)

// computes dotproduct of two vectors, which are given as points 
// struct point vector1, vector2, dest, 
// result goes into double dest
#define DOTPRODUCTp(vector1, vector2, dest)\
  DOTPRODUCT( (vector1).x, (vector1).y, (vector1).z,\
        (vector2).x, (vector2).y, (vector2).z,\
        (dest) )

#define DOTPRODUCT(x1, y1, z1, x2, y2, z2, dest)\
\
  dest = (x1) * (x2) + (y1) * (y2) + (z1) * (z2);

// normalizes vector dest
// struct point dest
// result returned in dest
// must declare a double variable called 'length' somewhere inside the scope of the NORMALIZE macrp
// macro will change that variable
#define pNORMALIZE(dest)\
\
  length = sqrt((dest).x * (dest).x + (dest).y * (dest).y + (dest).z * (dest).z);\
  (dest).x /= length;\
  (dest).y /= length;\
  (dest).z /= length;

// copies vector source to vector dest
// struct point source,dest
#define pCPY(source,dest)\
\
  (dest).x = (source).x;\
  (dest).y = (source).y;\
  (dest).z = (source).z;
  
// assigns values x,y,z to point vector dest
// struct point dest
// double x,y,z
#define pMAKE(x,y,z,dest)\
\
  (dest).(x) = (x);\
  (dest).(y) = (y);\
  (dest).(z) = (z);

// sums points src1 and src2 to dest
// struct point src1,src2,dest
#define pSUM(src1,src2,dest)\
\
  (dest).x = (src1).x + (src2).x;\
  (dest).y = (src1).y + (src2).y;\
  (dest).z = (src1).z + (src2).z;

// dest = src2 - src1
// struct point src1,src2,dest
#define pDIFFERENCE(src1,src2,dest)\
\
  (dest).x = (src1).x - (src2).x;\
  (dest).y = (src1).y - (src2).y;\
  (dest).z = (src1).z - (src2).z;

// mulitplies components of point src by scalar and returns the result in dest
// struct point src,dest
// double scalar
#define pMULTIPLY(src,scalar,dest)\
\
  (dest).x = (src).x * (scalar);\
  (dest).y = (src).y * (scalar);\
  (dest).z = (src).z * (scalar);

#endif


class MassSpring {
  public:  
  const double restLength = 1.0 / 7.0;
  const double restLengthDiagonal2D = restLength * sqrt(2.0);
  const double restLengthDiagonal3D = sqrt(restLength * restLength + restLengthDiagonal2D * restLengthDiagonal2D);
  const double restLengthBend = 2.0 / 7.0;

  // Calculate hook force between point a and its neighbour point b
  void calculateHookForce(const point& a, const point& b, double kHook, point& tempForce)
  {
    // Initialize temp force
    tempForce.x = 0;
    tempForce.y = 0;
    tempForce.z = 0;
    point vectBA;
    point normalizedVectBA;
    double vectBALength;
    double elasticForce;
    // Calculate vector pointing from b to a
    pDIFFERENCE(a, b, vectBA);
    vectBALength = sqrt(vectBA.x * vectBA.x + vectBA.y * vectBA.y + vectBA.z * vectBA.z);
    // Normalized vector vectBA
    normalizedVectBA.x = vectBA.x / vectBALength;
    normalizedVectBA.y = vectBA.y / vectBALength;
    normalizedVectBA.z = vectBA.z / vectBALength;
    // Elastic force exerted on a
    elasticForce = (-kHook) * (vectBALength - restLength);
    pMULTIPLY(normalizedVectBA, elasticForce, tempForce);
  }

  // Calculate hook force between point a and its diagonal neighbour point b on surface
  void calculateHookForceDiagonal2D(const point& a, const point& b, double kHook, point& tempForce)
  {
    // Initialize temp force
    tempForce.x = 0;
    tempForce.y = 0;
    tempForce.z = 0;
    point vectBA;
    point normalizedVectBA;
    double vectBALength;
    double elasticForce;
    // Calculate vector pointing from b to a
    pDIFFERENCE(a, b, vectBA);
    vectBALength = sqrt(vectBA.x * vectBA.x + vectBA.y * vectBA.y + vectBA.z * vectBA.z);
    // Normalized vector vectBA
    normalizedVectBA.x = vectBA.x / vectBALength;
    normalizedVectBA.y = vectBA.y / vectBALength;
    normalizedVectBA.z = vectBA.z / vectBALength;
    // Elastic force exerted on a
    elasticForce = (-kHook) * (vectBALength - restLengthDiagonal2D);
    pMULTIPLY(normalizedVectBA, elasticForce, tempForce);
  }

  // calculate hook force between point a and its diagonal neighbour point b in cube
  void calculateHookForceDiagonal3D(const point& a, const point& b, double kHook, point& tempForce)
  {
    // Initialize temp force
    tempForce.x = 0;
    tempForce.y = 0;
    tempForce.z = 0;
    point vectBA;
    point normalizedVectBA;
    double vectBALength;
    double elasticForce;
    // Calculate vector pointing from b to a
    pDIFFERENCE(a, b, vectBA);
    vectBALength = sqrt(vectBA.x * vectBA.x + vectBA.y * vectBA.y + vectBA.z * vectBA.z);
    // Normalized vector vectBA
    normalizedVectBA.x = vectBA.x / vectBALength;
    normalizedVectBA.y = vectBA.y / vectBALength;
    normalizedVectBA.z = vectBA.z / vectBALength;
    // Elastic force exerted on a
    elasticForce = (-kHook) * (vectBALength - restLengthDiagonal3D);
    pMULTIPLY(normalizedVectBA, elasticForce, tempForce);
  }

  // Calculate hook force between point a and its the next neighbour point b
  void calculateHookForceBend(const point& a, const point& b, double kHook, point& tempForce)
  {
    // Initialize temp force
    tempForce.x = 0;
    tempForce.y = 0;
    tempForce.z = 0;
    point vectBA;
    point normalizedVectBA;
    double vectBALength;
    double elasticForce;
    // Calculate vector pointing from b to a
    pDIFFERENCE(a, b, vectBA);
    vectBALength = sqrt(vectBA.x * vectBA.x + vectBA.y * vectBA.y + vectBA.z * vectBA.z);
    // Normalized vector vectBA
    normalizedVectBA.x = vectBA.x / vectBALength;
    normalizedVectBA.y = vectBA.y / vectBALength;
    normalizedVectBA.z = vectBA.z / vectBALength;
    // Elastic force exerted on a
    elasticForce = (-kHook) * (vectBALength - restLengthBend);
    pMULTIPLY(normalizedVectBA, elasticForce, tempForce);
  }

  // Calculate hook force for collision springs
  void calculateHookForceCollision(const point& a, const point& b, double kHook, point& n, point& tempForce)
  {
    // Initialize temp force
    tempForce.x = 0;
    tempForce.y = 0;
    tempForce.z = 0;
    point vectBA;
    double elasticForce;
    // Calculate vector pointing from b to a
    pDIFFERENCE(a, b, vectBA);
    pMULTIPLY(vectBA, -kHook, vectBA);
    DOTPRODUCTp(vectBA, n, elasticForce);
    pMULTIPLY(n, elasticForce, tempForce);
  }

  // Calculate damping force between point a and point b
  void calculateDampingForce(const point& a, const point& b, const point& va, const point& vb, double kDamp, point& tempForce)
  {
    // Initialize temp force
    tempForce.x = 0;
    tempForce.y = 0;
    tempForce.z = 0;
    point vectBA;
    point vectVelocity;
    point normalizedVectBA;
    double vectBALength;
    double elasticForce;
    // Calculate vector pointing from b to a
    pDIFFERENCE(a, b, vectBA);
    vectBALength = sqrt(vectBA.x * vectBA.x + vectBA.y * vectBA.y + vectBA.z * vectBA.z);
    // va and vb are velocities of points a and b
    pDIFFERENCE(va, vb, vectVelocity);
    elasticForce = (-kDamp) * ((vectVelocity.x * vectBA.x + vectVelocity.y * vectBA.y + vectVelocity.z * vectBA.z) / vectBALength);
    // Normalized vector vectBA
    normalizedVectBA.x = vectBA.x / vectBALength;
    normalizedVectBA.y = vectBA.y / vectBALength;
    normalizedVectBA.z = vectBA.z / vectBALength;
    pMULTIPLY(normalizedVectBA, elasticForce, tempForce);
  }

  // Calculate the damping force for collision springs
  void calculateDampingForceCollision(const point& a, const point& b, const point& va, const point& vb, double kDamp, point& n, point& tempForce)
  {
    // Initialize temp force
    tempForce.x = 0;
    tempForce.y = 0;
    tempForce.z = 0;
    point vectBA;
    point vectVelocity;
    point normalizedVectBA;
    double vectBALength;
    double elasticForce;
    double cosF;
    // Calculate vector pointing from b to a
    pDIFFERENCE(a, b, vectBA);
    vectBALength = sqrt(vectBA.x * vectBA.x + vectBA.y * vectBA.y + vectBA.z * vectBA.z);
    // va and vb are velocities of points a and b
    pDIFFERENCE(va, vb, vectVelocity);
    elasticForce = (-kDamp) * ((vectVelocity.x * vectBA.x + vectVelocity.y * vectBA.y + vectVelocity.z * vectBA.z) / vectBALength);
    // Normalized vector vectBA
    normalizedVectBA.x = vectBA.x / vectBALength;
    normalizedVectBA.y = vectBA.y / vectBALength;
    normalizedVectBA.z = vectBA.z / vectBALength;
    pMULTIPLY(normalizedVectBA, elasticForce, tempForce);
    DOTPRODUCTp(tempForce, n, cosF);
    pMULTIPLY(n, cosF, tempForce);
  }

  // Calculate structual force between point p at (i, j, k) and its neighbour point b
  void calculateStructualForce(world *jello, int i, int j, int k, point& a)
  {
    point f;
    if (i > 0) // point p at (i, j, k) has its left neighbour (i-1, j, k), calculate the hook force and damping force from left
    {
      calculateHookForce(jello->p[i][j][k], jello->p[i - 1][j][k], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i - 1][j][k], jello->v[i][j][k], jello->v[i - 1][j][k], jello->dElastic, f);
      pSUM(f, a, a);
    }
    if (i < 7) // point p at (i, j, k) has its right neighbour (i+1, j, k), calculate the hook force and damping force from right
    {
      calculateHookForce(jello->p[i][j][k], jello->p[i + 1][j][k], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i + 1][j][k], jello->v[i][j][k], jello->v[i + 1][j][k], jello->dElastic, f);
      pSUM(f, a, a);
    }
    if (j > 0) // point p at (i, j, k) has its upper neighbour (i, j-1, k), calculate the hook force and damping force from upper
    {
      calculateHookForce(jello->p[i][j][k], jello->p[i][j - 1][k], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i][j - 1][k], jello->v[i][j][k], jello->v[i][j - 1][k], jello->dElastic, f);
      pSUM(f, a, a);
    }
    if (j < 7) // point p at (i, j, k) has its bottom neighbour (i, j+1, k), calculate the hook force and damping force from bottom
    {
      calculateHookForce(jello->p[i][j][k], jello->p[i][j + 1][k], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i][j + 1][k], jello->v[i][j][k], jello->v[i][j + 1][k], jello->dElastic, f);
      pSUM(f, a, a);
    }
    if (k > 0) // point p at (i, j, k) has its front neighbour (i, j, k-1), calculate the hook force and damping force from front
    {
      calculateHookForce(jello->p[i][j][k], jello->p[i][j][k-1], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i][j][k-1], jello->v[i][j][k], jello->v[i][j][k-1], jello->dElastic, f);
      pSUM(f, a, a);
    }
    if (k < 7) // point p at (i, j, k) has its behind neighbour (i, j, k+1), calculate the hook force and damping force from behind
    {
      calculateHookForce(jello->p[i][j][k], jello->p[i][j][k+1], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i][j][k+1], jello->v[i][j][k], jello->v[i - 1][j][k], jello->dElastic, f);
      pSUM(f, a, a);
    }
  }

  // Calculate shear force between point p at (i, j, k) and its diagonal neighbour point b
  void calculateShearForce(world *jello, int i, int j, int k, point& a)
  {
    point f;
    if (i > 0) 
    {
      // Diagonals on 2D surface
      if (k > 0) // point p at (i, j, k) has its neighbour (i-1, j, k-1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j][k - 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i - 1][j][k - 1], jello->v[i][j][k], jello->v[i - 1][j][k - 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (k < 7) // point p at (i, j, k) has its neighbour (i-1, j, k+1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j][k + 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i - 1][j][k + 1], jello->v[i][j][k], jello->v[i - 1][j][k + 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (j > 0) // point p at (i, j, k) has its neighbour (i-1, j-1, k), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j - 1][k], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k], jello->v[i][j][k], jello->v[i - 1][j - 1][k], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (j < 7) // point p at (i, j, k) has its neighbour (i-1, j+1, k), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j + 1][k], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k], jello->v[i][j][k], jello->v[i - 1][j + 1][k], jello->dElastic, f);
        pSUM(f, a, a);
      }
      // Diagonals in 3D cube
      if (j > 0 && k > 0) // point p at (i, j, k) has its neighbour (i-1, j-1, k-1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j - 1][k - 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k - 1], jello->v[i][j][k], jello->v[i - 1][j - 1][k - 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (j > 0 && k < 7) // point p at (i, j, k) has its neighbour (i-1, j-1, k+1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j - 1][k + 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k + 1], jello->v[i][j][k], jello->v[i - 1][j - 1][k + 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (j < 7 && k > 0) // point p at (i, j, k) has its neighbour (i-1, j+1, k-1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j + 1][k - 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k - 1], jello->v[i][j][k], jello->v[i - 1][j + 1][k - 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (j < 7 && k < 7) // point p at (i, j, k) has its neighbour (i-1, j+1, k+1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j + 1][k + 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k + 1], jello->v[i][j][k], jello->v[i - 1][j + 1][k + 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
    }
    if (i < 7)
    {
      // Diagonals on 2D surface
      if (k > 0) // point p at (i, j, k) has its neighbour (i+1, j, k-1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j][k - 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i + 1][j][k - 1], jello->v[i][j][k], jello->v[i + 1][j][k - 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (k < 7) // point p at (i, j, k) has its neighbour (i+1, j, k+1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j][k + 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i + 1][j][k + 1], jello->v[i][j][k], jello->v[i + 1][j][k + 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (j > 0) // point p at (i, j, k) has its neighbour (i+1, j-1, k), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j - 1][k], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k], jello->v[i][j][k], jello->v[i + 1][j - 1][k], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (j < 7) // point p at (i, j, k) has its neighbour (i+1, j+1, k), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j + 1][k], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k], jello->v[i][j][k], jello->v[i + 1][j + 1][k], jello->dElastic, f);
        pSUM(f, a, a);
      }
      // Diagonals in 3D cube
      if (j > 0 && k > 0) // point p at (i, j, k) has its neighbour (i+1, j-1, k-1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j - 1][k - 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k - 1], jello->v[i][j][k], jello->v[i + 1][j - 1][k - 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (j > 0 && k < 7) // point p at (i, j, k) has its neighbour (i+1, j-1, k+1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j - 1][k + 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k + 1], jello->v[i][j][k], jello->v[i + 1][j - 1][k + 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (j < 7 && k > 0) // point p at (i, j, k) has its neighbour (i+1, j+1, k-1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j + 1][k - 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k - 1], jello->v[i][j][k], jello->v[i + 1][j + 1][k - 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (j < 7 && k < 7) // point p at (i, j, k) has its neighbour (i+1, j+1, k+1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j + 1][k + 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k + 1], jello->v[i][j][k], jello->v[i + 1][j + 1][k + 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
    }
    if (j > 0)
    {
      // Diagonals on 2D surface
      if (k > 0) // point p at (i, j, k) has its neighbour (i, j-1, k-1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j - 1][k - 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i][j - 1][k - 1], jello->v[i][j][k], jello->v[i][j - 1][k - 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (k < 7) // point p at (i, j, k) has its neighbour (i, j-1, k+1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j - 1][k + 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i][j - 1][k + 1], jello->v[i][j][k], jello->v[i][j - 1][k + 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
    }
    if (j < 7) {
      // on 2D surface
      if (k > 0) // point p at (i, j, k) has its neighbour (i, j+1, k-1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j + 1][k - 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i][j + 1][k - 1], jello->v[i][j][k], jello->v[i][j + 1][k - 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
      if (k < 7) // point p at (i, j, k) has its neighbour (i, j+1, k+1), calculate the hook force and damping force
      {
        calculateHookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j + 1][k + 1], jello->kElastic, f);
        pSUM(f, a, a);
        calculateDampingForce(jello->p[i][j][k], jello->p[i][j + 1][k + 1], jello->v[i][j][k], jello->v[i][j + 1][k + 1], jello->dElastic, f);
        pSUM(f, a, a);
      }
    }
  }

  // Calculate bend force between point p at (i, j, k) and its following neighbour point b
  void calculateBendForce(world *jello, int i, int j, int k, point& a) {
    point f;
    if (i > 1) // point p at (i, j, k) has its following neighbour (i-2, j, k), calculate the hook force and damping force
    {
      calculateHookForceBend(jello->p[i][j][k], jello->p[i - 2][j][k], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i - 2][j][k], jello->v[i][j][k], jello->v[i - 2][j][k], jello->dElastic, f);
      pSUM(f, a, a);
    }
    if (i < 6) // point p at (i, j, k) has its following neighbour (i+2, j, k), calculate the hook force and damping force
    {
      calculateHookForceBend(jello->p[i][j][k], jello->p[i + 2][j][k], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i + 2][j][k], jello->v[i][j][k], jello->v[i + 2][j][k], jello->dElastic, f);
      pSUM(f, a, a);
    } 
    if (j > 1) // point p at (i, j, k) has its following neighbour (i, j-2, k), calculate the hook force and damping force
    {
      calculateHookForceBend(jello->p[i][j][k], jello->p[i][j - 2][k], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i][j - 2][k], jello->v[i][j][k], jello->v[i][j - 2][k], jello->dElastic, f);
      pSUM(f, a, a);
    }
    if (j < 6) // point p at (i, j, k) has its following neighbour (i, j+2, k), calculate the hook force and damping force
    {
      calculateHookForceBend(jello->p[i][j][k], jello->p[i][j + 2][k], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i][j + 2][k], jello->v[i][j][k], jello->v[i][j + 2][k], jello->dElastic, f);
      pSUM(f, a, a);
    }
    if (k > 1) // point p at (i, j, k) has its following neighbour (i, j, k-2), calculate the hook force and damping force
    {
      calculateHookForceBend(jello->p[i][j][k], jello->p[i][j][k - 2], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i][j][k - 2], jello->v[i][j][k], jello->v[i][j][k - 2], jello->dElastic, f);
      pSUM(f, a, a);
    }
    if (k < 6) // point p at (i, j, k) has its following neighbour (i, j, k+2), calculate the hook force and damping force
    {
      calculateHookForceBend(jello->p[i][j][k], jello->p[i][j][k + 2], jello->kElastic, f);
      pSUM(f, a, a);
      calculateDampingForce(jello->p[i][j][k], jello->p[i][j][k + 2], jello->v[i][j][k], jello->v[i][j][k + 2], jello->dElastic, f);
      pSUM(f, a, a);
    }
  }

  // Calculate external force field, add it to point p at (i, j, k) 
  void calculateExternalForce(world *jello, int x, int y, int z, point& a) {
    // External force index in resolution array
    int i, j, k;
    // Forces at 8 corners in a specific grid
    point f000, f001; 
    point f010, f011;
    point f100, f101; 
    point f110, f111;
    // External force position in grid
    double px, py, pz;
    // Force field grid
    double grid; 
    // External force field value
    point externalForce;
    externalForce.x = 0;
    externalForce.y = 0;
    externalForce.z = 0;

    i = int((jello->p[x][y][z].x + 2) * (jello->resolution - 1) / 4);
    j = int((jello->p[x][y][z].y + 2) * (jello->resolution - 1) / 4);
    k = int((jello->p[x][y][z].z + 2) * (jello->resolution - 1) / 4);

    // Check if the index is at the wall of the bounding box
    if (i == (jello->resolution - 1)) {
      i--;
    }
    if (j == (jello->resolution - 1)) {
      j--;
    }
    if (k == (jello->resolution - 1)) {
      k--;
    }
    // Check if the point is inside the bounding box, read the force field value
    if (((i >= 0) && (i <= jello->resolution - 1)) && ((j >= 0) && (j <= jello->resolution - 1)) && ((j >= 0) && (j <= jello->resolution - 1))) {
      f000 = jello->forceField[(i * jello->resolution * jello->resolution + j * jello->resolution + k)];
      f001 = jello->forceField[(i * jello->resolution * jello->resolution + j * jello->resolution + (k + 1))];
       
      f010 = jello->forceField[(i * jello->resolution * jello->resolution + (j + 1) * jello->resolution + k)];
      f011 = jello->forceField[(i * jello->resolution * jello->resolution + (j + 1) * jello->resolution + (k + 1))];
      
      f100 = jello->forceField[((i + 1) * jello->resolution * jello->resolution + j * jello->resolution + k)];
      f101 = jello->forceField[((i + 1) * jello->resolution * jello->resolution + j * jello->resolution + (k + 1))];
      
      f110 = jello->forceField[((i + 1) * jello->resolution * jello->resolution + (j + 1) * jello->resolution + k)];
      f111 = jello->forceField[((i + 1) * jello->resolution * jello->resolution + (j + 1) * jello->resolution + (k + 1))];

      // 3D interpolation
      grid = 1.0 * 4 / (jello->resolution - 1);
      px = (jello->p[x][y][z].x - (-2 + 1.0 * 4 * i / (jello->resolution - 1))) / grid;
      py = (jello->p[x][y][z].y - (-2 + 1.0 * 4 * j / (jello->resolution - 1))) / grid;
      pz = (jello->p[x][y][z].z - (-2 + 1.0 * 4 * k / (jello->resolution - 1))) / grid;

      pMULTIPLY(f000, (1 - px) * (1 - py) * (1 - pz), f000);
      pMULTIPLY(f001, (1 - px) * (1 - py) * pz, f001);
      pMULTIPLY(f010, (1 - px) * py * (1 - pz), f010);
      pMULTIPLY(f011, (1 - px) * py * pz, f011);
      pMULTIPLY(f100, px * (1 - py) * (1 - pz), f100);
      pMULTIPLY(f101, px * (1 - py) * pz, f101);
      pMULTIPLY(f110, px * py * (1 - pz), f110);
      pMULTIPLY(f111, px * py * pz, f111);

      pSUM(externalForce, f000, externalForce);
      pSUM(externalForce, f001, externalForce);
      pSUM(externalForce, f010, externalForce);
      pSUM(externalForce, f011, externalForce);
      pSUM(externalForce, f100, externalForce);
      pSUM(externalForce, f101, externalForce);
      pSUM(externalForce, f110, externalForce);
      pSUM(externalForce, f111, externalForce);
      a.x = externalForce.x;
      a.y = externalForce.y;
      a.z = externalForce.z;
    }
  }

  // Calculate collision response force by the bounding box walls
  void calculateCollisionForce(world *jello, int i, int j, int k, point& a) {
    // Collision detection, check if point (i, j, k) is outside the bounding box
    if ((jello->p[i][j][k].x <= -2) || (jello->p[i][j][k].x >= 2) || 
      (jello->p[i][j][k].y <= -2) || (jello->p[i][j][k].y >= 2) || 
      (jello->p[i][j][k].z <= -2) || (jello->p[i][j][k].z >= 2)) {
      point f;
      // Initialize the normal vector of the wall
      point normal;
      normal.x = 0;
      normal.y = 0;
      normal.z = 0;
      // Initialize the collision point on the wall
      point wall;
      wall = jello->p[i][j][k];

      point vWall;
      vWall.x = 0;
      vWall.y = 0;
      vWall.z = 0;
      // Calculate the force direction and the collision point on the wall using penalty method
      if(jello->p[i][j][k].x <= -2)
      {
        normal.x = 1;
        wall.x = -2;
      }
      if(jello->p[i][j][k].x >= 2)
      {
        normal.x = -1;
        wall.x = 2;
      }
      if(jello->p[i][j][k].y <= -2)
      {
        normal.y = 1;
        wall.y = -2;
      }
      if(jello->p[i][j][k].y >= 2)
      {
        normal.y = -1;
        wall.y = 2;
      } 
      if(jello->p[i][j][k].z <= -2)
      {
        normal.z = 1;
        wall.z = -2;
      }
      if(jello->p[i][j][k].z >= 2)
      {
        normal.z = -1;
        wall.z = 2;
      }
      calculateHookForceCollision(jello->p[i][j][k], wall, jello->kCollision, normal, f);
      pSUM(f, a, a);
      calculateDampingForceCollision(jello->p[i][j][k], wall, jello->v[i][j][k], vWall, jello->dCollision, normal, f);
      pSUM(f, a, a);
    }
  }

  // Check the side of the plane and the given point
  void inclinedPlaneSide(world *jello, point& p, bool& sameSide, point& normal) {
    double point; 
    double plane; 
    int pointSide;
    int planeSide;
    
    // Inclined plane is presented in the bounding box
    if (jello->incPlanePresent == 1)  {
      // Check the side of the given point
      point = jello->a * p.x + jello->b * p.y + jello->c * p.z + jello->d;
      if (point > 0) {
        pointSide = 1;
      }
      else if (point < 0) {
        pointSide = -1;
      }
      
      // Check the orientation of the plane normal
      normal.x = jello->a;
      normal.y = jello->b;
      normal.z = jello->c;
      plane = jello->a * normal.x + jello->b * normal.y + jello->c * normal.z + jello->d;
      if (plane > 0) {
        planeSide = 1;
      }
      else if (plane < 0) {
        planeSide = -1;
      }

      // the plane normal and the point is on the same side
      if ((planeSide * pointSide)  == 1) {
        sameSide = true;
      }
      // the plane normal and the point is on different sides
      if ((planeSide * pointSide) == -1) {
        sameSide = false;
      }
    }
  }

  void inclinedPlaneCollisionForce(world * jello, int i, int j, int k, const bool& s, const point& pn, point &a) {
    // the point collided with the plane
    if (s == false) {
      double distance; 
      point normal;
      double normalLength;
      point normalizedN;
      point pointOnPlane;
      point f;
      point vPlane;
      vPlane.x = 0;
      vPlane.y = 0;
      vPlane.z = 0;
      // D = (ax + by + cz + d) / sqrt(a*a + b*b + c*c)
      distance = (jello->a * jello->p[i][j][k].x + jello->b * jello->p[i][j][k].y + jello->c * jello->p[i][j][k].z + jello->d);
      distance = distance / sqrt(jello->a * jello->a + jello->b * jello->b + jello->c * jello->c);
      
      normalLength = sqrt(pn.x * pn.x + pn.y * pn.y + pn.z * pn.z);
      // Normalized normal, get the orientation of the normal
      normalizedN.x = (pn.x / normalLength);
      normalizedN.y = (pn.y / normalLength);
      normalizedN.z = (pn.z / normalLength);
      pMULTIPLY(normalizedN, distance, normalizedN);
      // Calculate the point on plane
      pDIFFERENCE(jello->p[i][j][k], normalizedN, pointOnPlane);
      // Calculate the direction
      normal.x = (pn.x / normalLength);
      normal.y = (pn.y / normalLength);
      normal.z = (pn.z / normalLength);
      calculateHookForceCollision(jello->p[i][j][k], pointOnPlane, jello->kCollision, normal, f);
      pSUM(f, a, a);
      calculateDampingForceCollision(jello->p[i][j][k], pointOnPlane, jello->v[i][j][k], vPlane, jello->dCollision, normal, f);
      pSUM(f, a, a);
    }
  }
};
