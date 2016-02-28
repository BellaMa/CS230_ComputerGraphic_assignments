 /**
只修改了vertices有关的那俩vector
 */
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <stack>
#include <algorithm>
#include "math.h"
#include "minigl.h"
using namespace std;

const MGLfloat Z_MAX=2;


int current_mode = -1;
MGLpoly_mode poly_mode;
MGLpixel c_color;

class Vertex_3D
{
public:
  MGLfloat x, y, z, w;
  MGLpixel color;
  
  Vertex_3D()
  :x(0), y(0), z(0), w(1), color(c_color)
  {}
  Vertex_3D(MGLfloat _x, MGLfloat _y, MGLfloat _z)
  :x(_x), y(_y), z(_z), w(1), color(c_color)
  {}
  Vertex_3D(MGLfloat _x, MGLfloat _y, MGLfloat _z, MGLpixel _color)
  :x(_x), y(_y), z(_z), w(1), color(_color)
  {}
  Vertex_3D(MGLfloat * v)
  :color(c_color)
  {          
    x = v[0];
    y = v[1];
    z = v[2];
    w = v[3];         
  }
  Vertex_3D& operator= (const Vertex_3D & rhs){
    if (this != &rhs){
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        w= rhs.w;
        color = rhs.color;
    }
    return *this;
  }        

};

class Matrix{
public:
  MGLfloat t[16];
  
  Matrix() {
            t[0] = 1; t[4] = 0; t[8] = 0 ; t[12] = 0;
            t[1] = 0; t[5] = 1; t[9] = 0 ; t[12] = 0;
            t[2] = 0; t[6] = 0; t[10] = 1; t[14] = 0;
            t[3] = 0; t[7] = 0; t[11] = 0; t[15] = 1;
  }

  Matrix( MGLfloat * matrix) {
    for(int i=0;i<16;i++) {
        t[i]=matrix[i];
    }
  }
  Matrix(Matrix * matrix){
    for(int i=0;i<16;i++) {
        t[i]=matrix->t[i];
    }
  }
  Matrix(MGLfloat _t0, MGLfloat _t1, MGLfloat _t2, MGLfloat _t3,MGLfloat _t4,MGLfloat _t5,MGLfloat _t6,MGLfloat _t7,
               MGLfloat _t8,MGLfloat _t9,MGLfloat _t10,MGLfloat _t11,MGLfloat _t12,MGLfloat _t13,MGLfloat _t14,MGLfloat _t15) 
  {
    t[0]=_t0, t[1]=_t1,t[2]=_t2,t[3]=_t3,t[4]=_t4,t[5]=_t5,t[6]=_t6,t[7]=_t7,t[8]=_t8,t[9]=_t9,t[10]=_t10,t[11]=_t11,t[12]=_t12,t[13]=_t13,t[14]=_t14,t[15]=_t15;
  }

  Matrix operator * (Matrix& matrix){
    Matrix a;
    int z = 0;
    for(int i=0; i<=12;i=i+4) {
        for(int j=0; j<=3;j++) {
            a.t[z]=this->t[j] * matrix.t[i] + this->t[j+4] * matrix.t[i+1] + this->t[j+8] * matrix.t[i+2] + this->t[j+12] * matrix.t[i+3];
            z++;
        }
    }
    return a;
  };
  
  Vertex_3D operator * (Vertex_3D& rhs){
    MGLfloat temp[0];
    for(int i = 0 ; i < 4; ++i)
    {
      temp[i] = t[i]*rhs.x + t[i+4]*rhs.y + t[i+8]*rhs.z + t[i+12]*rhs.w;
    }
    MGLfloat w1 = temp[3];
    for (int i = 0; i < 3; ++i)
        temp[i] = temp[i] / w1;
    temp[3] = 1;
    Vertex_3D v(temp);   
    v.color = rhs.color;
    return v;
  }
};

class Triangle{
    public:
        Vertex_3D tri_a, tri_b, tri_c;
};

vector<Vertex_3D >  obj_vertices;
Matrix * currentMatrix_pointer;
vector<Triangle> vertices_tri;
stack<Matrix> modelviewStack, projectionStack;

void screenCoordinates(Vertex_3D& v, MGLsize w, MGLsize h){
    v.x = (v.x / 2.0 + 0.5) * (w + 1) + 0.5;
    v.y = (v.y / 2.0 + 0.5) * (h + 1) + 0.5;
}

MGLfloat lineEquation(Vertex_3D &v1, Vertex_3D &v2, Vertex_3D &v3)
{
    MGLfloat dy1 = v1.y - v2.y;
    MGLfloat dx1 = v3.x - v1.x;
    MGLfloat dx2 = v2.x - v1.x;
    MGLfloat dy2 = v3.y - v1.y;
    return (dy1 * dx1 + dx2 * dy2);
}
/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(char* description) {
    printf("%s\n", description);
    exit(1);
}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */

void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
    vector<MGLfloat> z_buffer;
    for(int i = 0; i < width * height; ++i){
        z_buffer.push_back(Z_MAX);
    } 

    Vertex_3D curr, a_tri,b_tri,c_tri;
    MGLfloat alpha,beta,gamma;
    

    for(int i = 0; i < vertices_tri.size(); ++i)
    {
        a_tri = vertices_tri.at(i).tri_a;
        b_tri = vertices_tri.at(i).tri_b;
        c_tri = vertices_tri.at(i).tri_c;

        screenCoordinates(a_tri, width, height);
        screenCoordinates(b_tri, width, height);
        screenCoordinates(c_tri, width, height);

        for(int x = 0; x < width; ++x){
            for(int y = 0; y < height; ++y){
                curr.x = (MGLfloat) x;
                curr.y = (MGLfloat) y;


                MGLfloat bc_curr = lineEquation(b_tri,c_tri,curr);
                MGLfloat bc =  lineEquation(b_tri,c_tri,a_tri);
                alpha = bc_curr / bc;

                MGLfloat ca_curr = lineEquation(c_tri,a_tri,curr);
                MGLfloat ca = lineEquation(c_tri,a_tri,b_tri);
                beta = ca_curr / ca;

                MGLfloat ab_curr = lineEquation(a_tri,b_tri,curr);
                MGLfloat ab = lineEquation(a_tri,b_tri,c_tri);
                gamma = ab_curr / ab;

                if( alpha >0 && beta >= 0 && gamma >= 0){
                    int p = y * width + x;
                    curr.z = alpha * a_tri.z + beta * b_tri.z + gamma * c_tri.z; 
                    if(curr.z < z_buffer.at(p)){
                        z_buffer.at(p) = curr.z;

                        MGLbyte color_red = alpha * MGL_GET_RED(a_tri.color) + beta * MGL_GET_RED(b_tri.color) + gamma * MGL_GET_RED(c_tri.color);
                        MGLbyte color_green = alpha * MGL_GET_GREEN(a_tri.color) + beta * MGL_GET_GREEN(b_tri.color) + gamma * MGL_GET_GREEN(c_tri.color);
                        MGLbyte color_blue = alpha * MGL_GET_BLUE(a_tri.color) + beta * MGL_GET_BLUE(b_tri.color) + gamma * MGL_GET_BLUE(c_tri.color);
                        
                        mglColor(color_red,color_green,color_blue, MGLpixel data[p]);

                    }

                }
            }
        }

    }

}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
  poly_mode = mode;
}

/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
  int i;
  switch (poly_mode) {
    case MGL_TRIANGLES:
      for(i = 0; i < obj_vertices.size(); i =i + 3){
        Triangle temp;
        temp.tri_a =  obj_vertices[i];
        temp.tri_b =  obj_vertices[i + 1];
        temp.tri_c =  obj_vertices[i + 2];
        vertices_tri.push_back(temp);
      }
      obj_vertices.clear();
      break;
      
    case MGL_QUADS:
      for(i = 0; i < obj_vertices.size(); i =i + 4){
            Triangle temp;
            temp.tri_a =  obj_vertices[i];
            temp.tri_b =  obj_vertices[i + 1];
            temp.tri_c =  obj_vertices[i + 2];
            vertices_tri.push_back(temp);

            temp.tri_a =  obj_vertices[i];
            temp.tri_b =  obj_vertices[i + 2];
            temp.tri_c =  obj_vertices[i + 3];
            vertices_tri.push_back(temp);
      }
      obj_vertices.clear();
      break;
  }
}

/**
 * Specify a two-dimensional Vertex_3D; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
  mglVertex3(x, y, 0);
}

/**
 * Specify a three-dimensional Vertex_3D.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
  Vertex_3D  v(x, y, z, c_color);
  v = projectionStack.top() * modelviewStack.top() * (v);
  obj_vertices.push_back(v);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
    Matrix temp;
    if(modelviewStack.empty()) modelviewStack.push(temp);
    if(projectionStack.empty()) projectionStack.push(temp);

    if(mode == MGL_MODELVIEW)
    {
        current_mode = mode;
        currentMatrix_pointer = & modelviewStack.top();
    }
    else if(mode == MGL_PROJECTION)
    {
        current_mode = mode;
        currentMatrix_pointer = & projectionStack.top();
    }
    else
        current_mode = -1; 
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
   Matrix temp(currentMatrix_pointer);

    if(current_mode == MGL_MODELVIEW){
      modelviewStack.push(temp);
      currentMatrix_pointer = & modelviewStack.top();
    } else if (current_mode == MGL_PROJECTION){
      projectionStack.push(temp);
      currentMatrix_pointer = & projectionStack.top();
    }

}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
    if(current_mode == MGL_MODELVIEW && !modelviewStack.empty()){        
        modelviewStack.pop();
        currentMatrix_pointer = & modelviewStack.top();
    }
    else if(current_mode == MGL_PROJECTION && !projectionStack.empty()){
        projectionStack.pop();
        currentMatrix_pointer = & projectionStack.top();       
    }
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix( MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix( Matrix & matrix)
{
  *currentMatrix_pointer = (*currentMatrix_pointer) * matrix;
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
  Matrix translate(1,0,0,0,0,1,0,0,0,0,1,0,x,y,z,1);
  mglMultMatrix(translate);
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
  MGLfloat t = sqrt(x*x + y*y + z*z);
  x = x / t;
  y = y / t;
  z = z / t;
  MGLfloat c = cos(angle * M_PI / 180.0),
           s = sin(angle * M_PI / 180.0);
  Matrix rotate ( x*x*(1-c)+c,y*x*(1-c)+z*s,z*x*(1-c)-y*s,0,x*y*(1-c)-z*s,y*y*(1-c)+c,z*y*(1-c)+x*s,0,x*z*(1-c)+y*s,y*z*(1-c)-x*s,z*z*(1-c)+c, 0,0,0,0,1);
  mglMultMatrix(rotate);
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
  Matrix scale(x,0,0,0,0,y,0,0,0,0,z,0,0,0,0,1);
  mglMultMatrix(scale);
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)

{ Matrix rustum(2.0*near/(right-left),0,0,0,0,2.0*near/(top-bottom),0,0,(right+left)/(right-left),(top + bottom)/(top-bottom),-(far+near)/(far-near),-1,0,0,-2.0*far*near/(far-near),0);
  mglMultMatrix(rustum);
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
  Matrix ortho(2.0/(right-left),0,0,0,0,2.0/(top-bottom),0,0,0,0,-2.0/(far-near),0,-(right+left)/(right-left),-(top+bottom)/(top-bottom),-(far+near)/(far-near),1);
  mglMultMatrix(ortho);
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLbyte red,
              MGLbyte green,
              MGLbyte blue,  MGLpixel *data)
{
  MGL_SET_BLUE(*data, blue);
  MGL_SET_RED(*data, red);
  MGL_SET_GREEN(*data, green);
}
