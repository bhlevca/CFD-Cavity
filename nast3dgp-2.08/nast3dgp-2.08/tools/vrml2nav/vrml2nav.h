#ifndef VRML2NAV_INCLUDED
#define VRML2NAV_INCLUDED
#include "tokenizer.h"

#define TOKEN_HEAD          1
#define TOKEN_COMMENT       2
#define TOKEN_ASCIITEXT     3
#define TOKEN_CONE          4
#define TOKEN_CUBE          5
#define TOKEN_CYLINDER      6
#define TOKEN_INDEXEDFACESET 7
#define TOKEN_INDEXEDLINESET 8
#define TOKEN_POINTSET      9
#define TOKEN_SPHERE        10
#define TOKEN_COORDINATE3   11
#define TOKEN_FONTSTYLE     12
#define TOKEN_INFO          13
#define TOKEN_LOD           14
#define TOKEN_MATERIAL      15
#define TOKEN_MATERIALBINDING 16
#define TOKEN_NORMAL        17
#define TOKEN_NORMALBINDING 18
#define TOKEN_TEXTURE2      19
#define TOKEN_TEXTURE2TRANSFORM 20
#define TOKEN_TEXTURECOORDINATE2 21
#define TOKEN_SHAPEHINTS    22
#define TOKEN_MATRIXTRANSFORM 23
#define TOKEN_ROTATION      24
#define TOKEN_SCALE         25
#define TOKEN_TRANSFORM     26
#define TOKEN_TRANSLATION   27
#define TOKEN_ORTHOGRAPHICCAMERA 28
#define TOKEN_PERSPECTIVECAMERA 29
#define TOKEN_DIRECTIONALLIGHT 30
#define TOKEN_POINTLIGHT    31
#define TOKEN_SPOTLIGHT     32
#define TOKEN_GROUP         33
#define TOKEN_SEPARATOR     34
#define TOKEN_SWITCH        35
#define TOKEN_TRANSFORMSEPARATOR 36
#define TOKEN_WWWANCHOR     37
#define TOKEN_WWWINLINE     38
#define TOKEN_STRING        39
#define TOKEN_SPACING       40
#define TOKEN_JUSTIFICATION 41
#define TOKEN_WIDTH         42
#define TOKEN_LCB           43
#define TOKEN_RCB           44
#define TOKEN_DAC           45
#define TOKEN_LEFT          46
#define TOKEN_RIGHT         47
#define TOKEN_CENTER        48
#define TOKEN_PARTS         49
#define TOKEN_BOTTOMRADIUS  50
#define TOKEN_HEIGHT        51
#define TOKEN_SIDES         52
#define TOKEN_BOTTOM        53
#define TOKEN_ALL           54
#define TOKEN_LEB           55
#define TOKEN_REB           56
#define TOKEN_NUMBER        57
#define TOKEN_COMMA         58
#define TOKEN_POINT         59
#define TOKEN_DEPTH         60
#define TOKEN_TOP           61
#define TOKEN_PIPE          62
#define TOKEN_RADIUS        63
#define TOKEN_ON            64
#define TOKEN_TRUE          65
#define TOKEN_FALSE         66
#define TOKEN_INTENSITY     67
#define TOKEN_COLOR         68 
#define TOKEN_DIRECTION     69
#define TOKEN_SIZE          70
#define TOKEN_FAMILY        71
#define TOKEN_SERIF         72
#define TOKEN_SANS          73
#define TOKEN_TYPEWRITER    74
#define TOKEN_STYLE         75
#define TOKEN_NONE          76
#define TOKEN_BOLD          77
#define TOKEN_ITALIC        78
#define TOKEN_COORDINDEX    79
#define TOKEN_MATERIALINDEX 80
#define TOKEN_NORMALINDEX   81
#define TOKEN_TEXTURECOORDINDEX 82
#define TOKEN_DEF           83
#define TOKEN_ROTATIONs     85
#define TOKEN_NAME          86
#define TOKEN_CENTERs       87
#define TOKEN_RANGE         88
#define TOKEN_AMBIENTCOLOR  89
#define TOKEN_DIFFUSECOLOR  90
#define TOKEN_SPECULARCOLOR 91
#define TOKEN_EMISSIVECOLOR 92
#define TOKEN_SHININESS     93
#define TOKEN_TRANSPARENCY  94
#define TOKEN_VALUE         95
#define TOKEN_DEFAULT       96
#define TOKEN_OVERALL       97
#define TOKEN_PER_PART      98
#define TOKEN_PER_PART_INDEXED 99
#define TOKEN_PER_FACE      100
#define TOKEN_PER_FACE_INDEXED 101
#define TOKEN_PER_VERTEX    102
#define TOKEN_PER_VERTEX_INDEXED 103
#define TOKEN_MATRIX        104
#define TOKEN_VECTOR        105
#define TOKEN_POSITION      106
#define TOKEN_ORIENTATION   107
#define TOKEN_FOCALDISTANCE 108
#define TOKEN_HEIGHTANGLE   109
#define TOKEN_LOCATION      112
#define TOKEN_STARTINDEX    113
#define TOKEN_NUMPOINTS     114
#define TOKEN_SCALEFACTOR   115
#define TOKEN_VERTEXORDERING 116
#define TOKEN_SHAPETYPE     117
#define TOKEN_FACETYPE      118
#define TOKEN_CREASEANGLE   119
#define TOKEN_UNKNOWN_ORDERING 120
#define TOKEN_CLOCKWISE     121
#define TOKEN_COUNTERCLOCKWISE 122
#define TOKEN_UNKNOWN_SHAPE_TYPE 123
#define TOKEN_SOLID         124
#define TOKEN_UNKNOWN_FACE_TYPE 125
#define TOKEN_CONVEX        126





class Vector3 {
public:
    double x,y,z;
    void operator=(Vector3& V) {x=V.x;y=V.y;z=V.z;}
    void matmul(double mat[3][3]) {
        double nx,ny,nz;
        nx=mat[0][0]*x+mat[0][1]*y+mat[0][2]*z;
        ny=mat[1][0]*x+mat[1][1]*y+mat[1][2]*z;
        nz=mat[2][0]*x+mat[2][1]*y+mat[2][2]*z;
        x=nx;y=ny;z=nz;
    }
};


#define JUSTIFICATION_LEFT 0
#define JUSTIFICATION_RIGHT 1
#define JUSTIFICATION_CENTER 2
class CAsciiText {
public:
    char*   string;
    double  spacing;
    int     justification;
    double  width;
            CAsciiText() {string=new char[1];string[0]=0;spacing=1;justification=JUSTIFICATION_LEFT;width=0;}
           ~CAsciiText() {if(string) delete string;}
};

#define PARTS_CONE_NONE 0
#define PARTS_CONE_SIDES 1
#define PARTS_CONE_BOTTOM 2
#define PARTS_CONE_ALL (1|2)
class CCone {
public:
    int     parts;
    double  bottomRadius;
    double  height;
            CCone() {parts=PARTS_CONE_ALL;bottomRadius=1;height=2;}
};

class CCoordinate3 {
public:
    int     size;
    Vector3* data;
            CCoordinate3() {size=0;data=0;}
           ~CCoordinate3() {if(data) delete data;}
    void    Dump() {
                for(int i=0;i<size;i++) printf("%f %f %f\n",data[i].x,data[i].y,data[i].z);
                }
};

class CCube {
public:
    double  width,height,depth;
            CCube() {width=height=depth=2;}
};


#define PARTS_CYLINDER_NONE 0
#define PARTS_CYLINDER_SIDES 1
#define PARTS_CYLINDER_TOP 2
#define PARTS_CYLINDER_BOTTOM 4
#define PARTS_CYLINDER_ALL (1|2|4)
class CCylinder {
public:
    int     parts;
    double  radius;
    double  height;
            CCylinder() {parts=PARTS_CYLINDER_ALL;radius=1;height=2;}
};

class CDirectionalLight {
public:
    int     on;
    double  intensity;
    double  color[3];
    Vector3 direction;
            CDirectionalLight() {on=1;intensity=1;color[0]=color[1]=color[2]=1;direction.x=direction.y=0;direction.z=-1;}
};

#define FAMILY_SERIF 0
#define FAMILY_SANS 1
#define FAMILY_TYPEWRITER 2
#define STYLE_NONE 0
#define STYLE_BOLD 1
#define STYLE_ITALIC 2
class CFontStyle {
public:
    double  size;
    int     family;
    int     style;
            CFontStyle() {size=10;family=FAMILY_SERIF;style=STYLE_NONE;}
};

class CIndexedFaceSet {
public:
    int     coordIndex_size;
    int     materialIndex_size;
    int     normalIndex_size;
    int     textureCoordIndex_size;
    long*   coordIndex;
    long*   materialIndex;
    long*   normalIndex;
    long*   textureCoordIndex;
            CIndexedFaceSet() {
                coordIndex_size=materialIndex_size=normalIndex_size=textureCoordIndex_size=1;
                coordIndex=new long[1];coordIndex[0]=0;
                materialIndex=new long[1];materialIndex[0]=-1;
                normalIndex=new long[1];normalIndex[0]=-1;
                textureCoordIndex=new long[1];textureCoordIndex[0]=-1;
                }
           ~CIndexedFaceSet() {
                if(coordIndex) delete coordIndex;
                if(materialIndex) delete materialIndex;
                if(normalIndex) delete normalIndex;
                if(textureCoordIndex) delete textureCoordIndex;
            }
    void    CloseCoordIndex();
};

class CIndexedLineSet {
public:
    int     coordIndex_size;
    int     materialIndex_size;
    int     normalIndex_size;
    int     textureCoordIndex_size;
    long*   coordIndex;
    long*   materialIndex;
    long*   normalIndex;
    long*   textureCoordIndex;
            CIndexedLineSet() {
                coordIndex_size=materialIndex_size=normalIndex_size=textureCoordIndex_size=1;
                coordIndex=new long[1];coordIndex[0]=0;
                materialIndex=new long[1];materialIndex[0]=-1;
                normalIndex=new long[1];normalIndex[0]=-1;
                textureCoordIndex=new long[1];textureCoordIndex[0]=-1;
                }
           ~CIndexedLineSet() {
                if(coordIndex) delete coordIndex;
                if(materialIndex) delete materialIndex;
                if(normalIndex) delete normalIndex;
                if(textureCoordIndex) delete textureCoordIndex;
            }
};

class CInfo {
public:
    char*   string;
            CInfo() {string=0;}
           ~CInfo() {if(string) delete string;}
};

class CLOD {
public:
    int     range_count;
    double*  range;
    double   center[3];
            CLOD() {range=0;center[0]=center[1]=center[2]=0;}
           ~CLOD() {if(range) delete range;}
};

class CMaterial {
public:    
    double   ambientColor[3];
    double   diffuseColor[3];
    double   specularColor[3];
    double   emissiveColor[3];
    double   shininess;
    double   transparency;
            CMaterial() {
                ambientColor[0]=ambientColor[1]=ambientColor[2]=0.2;
                diffuseColor[0]=diffuseColor[1]=diffuseColor[2]=0.8;
                specularColor[0]=specularColor[1]=specularColor[2]=0.0;
                emissiveColor[0]=emissiveColor[1]=emissiveColor[2]=0.0;
                shininess=0.2;
                transparency=0.0;
            }
};


#define VALUE_DEFAULT 1
#define VALUE_OVERALL 2
#define VALUE_PER_PART 3
#define VALUE_PER_PART_INDEXED 4
#define VALUE_PER_FACE 5
#define VALUE_PER_FACE_INDEXED 6
#define VALUE_PER_VERTEX 7
#define VALUE_PER_VERTEX_INDEXED 8
class CMaterialBinding {
public:
    int     value;
            CMaterialBinding() {value=VALUE_DEFAULT;}
};

class CMatrixTransform {
public:
    double   matrix[4][4];
            CMatrixTransform() {for(int i=0;i<4;i++) for(int j=0;j<4;j++) matrix[i][j]=(i==j?1:0);}
};

class CNormal {
public: 
  int     size;
  Vector3* data;
  double   vector[3];
  CNormal() {size=0; data=0; vector[0]=vector[1]=0;vector[2]=1;}
  ~CNormal() {if(data) delete data;}
  // CNormal() {vector[0]=vector[1]=0;vector[2]=1;}
};

class CNormalBinding {
public:
    int     value;
    CNormalBinding() {value=VALUE_DEFAULT;}
};

class COrthographicCamera {
public:
    double  position[3];
    double  orientation[4];
    double  focalDistance;
    double  height;
            COrthographicCamera() {
                position[0]=position[1]=orientation[0]=orientation[1]=orientation[3]=0;
                position[2]=orientation[2]=1;
                focalDistance=5;
                height=2;
            }
};

class CPerspectiveCamera {
public:
    double  position[3];
    double  orientation[4];
    double  focalDistance;
    double  heightAngle;
            CPerspectiveCamera() {
                position[0]=position[1]=orientation[0]=orientation[1]=orientation[3]=0;
                position[2]=orientation[2]=1;
                focalDistance=5;
                heightAngle=0.785398;
            }
};

class CPointLight {
public:
    int     on;
    double  intensity;
    double  color[3];
    double  location[3];
            CPointLight() {
                on=1;
                intensity=1;
                color[0]=color[1]=color[2]=1;
                location[0]=location[1]=0;location[2]=1;
            }
};

class CPointSet {
public:
    int     startIndex;
    int     numPoints;
            CPointSet() {startIndex=0;numPoints=-1;}
};

class CRotation {
public:
    Vector3 axis;
    double  angle;
    CRotation() {axis.x=axis.y=0;axis.z=1;angle=0;}
};

class CScale {
public:
    double  scaleFactor[3];
            CScale() {scaleFactor[0]=scaleFactor[1]=scaleFactor[2]=1;}
};


#define VERTEXORDERING_UNKNOWN_ORDERING 1
#define VERTEXORDERING_CLOCKWISE 2
#define VERTEXORDERING_COUNTERCLOCKWISE 3
#define SHAPETYPE_UNKNOWN_SHAPE_TYPE 1
#define SHAPETYPE_SOLID 2
#define FACETYPE_UNKNOWN_FACE_TYPE 1
#define FACETYPE_CONVEX 2
class CShapeHints {
public:
    int     vertexOrdering;
    int     shapeType;
    int     faceType;
    double  creaseAngle;
            CShapeHints() {
                vertexOrdering=VERTEXORDERING_UNKNOWN_ORDERING;
                shapeType=SHAPETYPE_UNKNOWN_SHAPE_TYPE;
                faceType=FACETYPE_CONVEX;
                creaseAngle=0.5;
            }
};

class ParseVRML {
public:
                ParseVRML(char* filename) {t=new Tokenizer(filename);}
                ParseVRML(FILE* f) {t=new Tokenizer(f);}
               ~ParseVRML() {delete t;}
    void        DoIt();
    void        DoItQuickAndDirty();

protected:
    Tokenizer*  t;


    int         PAsciiText(CAsciiText& C);
    int         PCone(CCone& C);    
    int         PCoordinate3(CCoordinate3& C);
    int         PCube(CCube& C);    
    int         PCylinder(CCylinder& C);
    int         PDirectionalLight(CDirectionalLight& C);
    int         PFontStyle(CFontStyle& C);
    int         PIndexedFaceSet(CIndexedFaceSet& C);
    int         PIndexedLineSet(CIndexedLineSet& C);
    int         PInfo(CInfo& C);
    int         PLOD(CLOD& C);
    int         PMaterial(CMaterial& C);
    int         PMaterialBinding(CMaterialBinding& C);
    int         PMatrixTransform(CMatrixTransform& C);
    int         PNormal(CNormal& C);
    int         PNormalBinding(CNormalBinding& C);
    int         POrthographicCamera(COrthographicCamera& C);
    int         PPerspectiveCamera(CPerspectiveCamera& C);
    int         PPointLight(CPointLight& C);
    int         PPointSet(CPointSet& C);
    int         PRotation(CRotation& C);
    int         PScale(CScale& C);
    int         PShapeHints(CShapeHints& C);

    int         PMFVec3f(int& count,Vector3** data);
    int         PMFLong(int& count,long** data);
    int         PMFFloat(int& count,double** data);
    int         ShowToken();
    int         GetToken();
    void        WriteHalfSpace(FILE* f,const Vector3& p1,const Vector3& p2,const Vector3& p3,const Vector3* pinside);

};







#endif
