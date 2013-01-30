#include "tokenizer.h"
#include "vrml2nav.h"
#define LISTNOFILE
#include "list.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include <vector>
#include <list>

#define sqr(x) ((x)*(x))
#define PI 3.1415926535897932384626433833


Token tok[]={
        {"#VRML V1.0 ascii",TOKEN_HEAD},
        {"#",TOKEN_COMMENT},
        {"AsciiText",TOKEN_ASCIITEXT},
        {"Cone",TOKEN_CONE},
        {"Cube",TOKEN_CUBE},
        {"Cylinder",TOKEN_CYLINDER},
        {"IndexedFaceSet",TOKEN_INDEXEDFACESET},
        {"IndexedLineSet",TOKEN_INDEXEDLINESET},
        {"PointSet",TOKEN_POINTSET},
        {"Sphere",TOKEN_SPHERE},
        {"Coordinate3",TOKEN_COORDINATE3},
        {"FontStyle",TOKEN_FONTSTYLE},
        {"Info",TOKEN_INFO},
        {"LOD",TOKEN_LOD},
        {"Material",TOKEN_MATERIAL},
        {"MaterialBinding",TOKEN_MATERIALBINDING},
        {"Normal",TOKEN_NORMAL},
        {"NormalBinding",TOKEN_NORMALBINDING},
        {"Texture2",TOKEN_TEXTURE2},
        {"Texture2Transform",TOKEN_TEXTURE2TRANSFORM},
        {"TextureCoordinate2",TOKEN_TEXTURECOORDINATE2},
        {"ShapeHints",TOKEN_SHAPEHINTS},
        {"MatrixTransform",TOKEN_MATRIXTRANSFORM},
        {"Rotation",TOKEN_ROTATION},
        {"Scale",TOKEN_SCALE},
        {"Transform",TOKEN_TRANSFORM},
        {"Translation",TOKEN_TRANSLATION},
        {"OrthographicCamera",TOKEN_ORTHOGRAPHICCAMERA},
        {"PerspectiveCamera",TOKEN_PERSPECTIVECAMERA},
        {"DirectionalLight",TOKEN_DIRECTIONALLIGHT},
        {"PointLight",TOKEN_POINTLIGHT},
        {"SpotLight",TOKEN_SPOTLIGHT},
        {"Group",TOKEN_GROUP},
        {"Separator",TOKEN_SEPARATOR},
        {"Switch",TOKEN_SWITCH},
        {"TransformSeparator",TOKEN_TRANSFORMSEPARATOR},
        {"WWWAnchor",TOKEN_WWWANCHOR},
        {"WWWInline",TOKEN_WWWINLINE},
        {"{",TOKEN_LCB},
        {"}",TOKEN_RCB},
        {"string",TOKEN_STRING},
        {"spacing",TOKEN_SPACING},
        {"justification",TOKEN_JUSTIFICATION},
        {"width",TOKEN_WIDTH},
        {"\"",TOKEN_DAC},
        {"LEFT",TOKEN_LEFT},
        {"RIGHT",TOKEN_RIGHT},
        {"CENTER",TOKEN_CENTER},
        {"parts",TOKEN_PARTS},
        {"bottomRadius",TOKEN_BOTTOMRADIUS},
        {"height",TOKEN_HEIGHT},
        {"SIDES",TOKEN_SIDES},
        {"BOTTOM",TOKEN_BOTTOM},
        {"ALL",TOKEN_ALL},
        {"[",TOKEN_LEB},
        {"]",TOKEN_REB},
        {"0",TOKEN_NUMBER},
        {"1",TOKEN_NUMBER},
        {"2",TOKEN_NUMBER},
        {"3",TOKEN_NUMBER},
        {"4",TOKEN_NUMBER},
        {"5",TOKEN_NUMBER},
        {"6",TOKEN_NUMBER},
        {"7",TOKEN_NUMBER},
        {"8",TOKEN_NUMBER},
        {"9",TOKEN_NUMBER},
        {"+",TOKEN_NUMBER},
        {"-",TOKEN_NUMBER},
        {",",TOKEN_COMMA},
        {"point",TOKEN_POINT},
        {"depth",TOKEN_DEPTH},
        {"TOP",TOKEN_TOP},
        {"|",TOKEN_PIPE},
        {"radius",TOKEN_RADIUS},
        {"on",TOKEN_ON},
        {"TRUE",TOKEN_TRUE},
        {"FALSE",TOKEN_FALSE},
        {"intensity",TOKEN_INTENSITY},
        {"color",TOKEN_COLOR},
        {"direction",TOKEN_DIRECTION},
        {"size",TOKEN_SIZE},
        {"family",TOKEN_FAMILY},
        {"SERIF",TOKEN_SERIF},
        {"SANS",TOKEN_SANS},
        {"TYPEWRITER",TOKEN_TYPEWRITER},
        {"style",TOKEN_STYLE},
        {"NONE",TOKEN_NONE},
        {"BOLD",TOKEN_BOLD},
        {"ITALIC",TOKEN_ITALIC},
        {"coordIndex",TOKEN_COORDINDEX},
        {"materialIndex",TOKEN_MATERIALINDEX},
        {"normalIndex",TOKEN_NORMALINDEX},
        {"textureCoordIndex",TOKEN_TEXTURECOORDINDEX},
        {"DEF",TOKEN_DEF},
        {"rotation",TOKEN_ROTATIONs},
        {"name",TOKEN_NAME},
    {"center",TOKEN_CENTERs},
    {"range",TOKEN_RANGE},
    {"ambientColor",TOKEN_AMBIENTCOLOR},
    {"diffuseColor",TOKEN_DIFFUSECOLOR},
    {"specularColor",TOKEN_SPECULARCOLOR},
    {"emissiveColor",TOKEN_EMISSIVECOLOR},
    {"shininess",TOKEN_SHININESS},
    {"transparency",TOKEN_TRANSPARENCY},
    {"value",TOKEN_VALUE},
    {"DEFAULT",TOKEN_DEFAULT},
    {"OVERALL",TOKEN_OVERALL},
    {"PER_PART",TOKEN_PER_PART},
    {"PER_PART_INDEXED",TOKEN_PER_PART_INDEXED},
    {"PER_FACE",TOKEN_PER_FACE},
    {"PER_FACE_INDEXED",TOKEN_PER_FACE_INDEXED},
    {"PER_VERTEX",TOKEN_PER_VERTEX},
    {"PER_VERTEX_INDEXED",TOKEN_PER_VERTEX_INDEXED},
    {"matrix",TOKEN_MATRIX},
    {"vector",TOKEN_VECTOR},
    {"position",TOKEN_POSITION},
    {"orientation",TOKEN_ORIENTATION},
    {"focalDistance",TOKEN_FOCALDISTANCE},
    {"heightAngle",TOKEN_HEIGHTANGLE},
    {"location",TOKEN_LOCATION},
    {"startIndex",TOKEN_STARTINDEX},
    {"numPoints",TOKEN_NUMPOINTS},
    {"scaleFactor",TOKEN_SCALEFACTOR},
    {"vertexOrdering",TOKEN_VERTEXORDERING},
    {"shapeType",TOKEN_SHAPETYPE},
    {"faceType",TOKEN_FACETYPE},
    {"creaseAngle",TOKEN_CREASEANGLE},
    {"UNKNOWN_ORDERING",TOKEN_UNKNOWN_ORDERING},
    {"CLOCKWISE",TOKEN_CLOCKWISE},
    {"COUNTERCLOCKWISE",TOKEN_COUNTERCLOCKWISE},
    {"UNKNOWN_SHAPE_TYPE",TOKEN_UNKNOWN_SHAPE_TYPE},
    {"SOLID",TOKEN_SOLID},
    {"UNKNOWN_FACE_TYPE",TOKEN_UNKNOWN_FACE_TYPE},
    {"CONVEX",TOKEN_CONVEX},

        {0,0}};




char* etchars="{[|";
char* syntax_error="Syntax error";


int ParseVRML::ShowToken() {
        int token;
        do {
                token=t->ShowToken(tok);
                if(token==TOKEN_COMMENT) {t->GetToken(tok);t->NewLine();}
        } while(token==TOKEN_COMMENT);
        return token;
}

int ParseVRML::GetToken() {
        int token;
        do {
                token=t->GetToken(tok);
                if(token==TOKEN_COMMENT) t->NewLine();
        } while(token==TOKEN_COMMENT);
        return token;
}


void ParseVRML::DoIt() {
        t->SetEndTokenChars(etchars);
        if(t->GetToken(tok)!=TOKEN_HEAD) t->Error("Not a VRML 1.0 file.");
        t->NewLine();
        CCone cl;
        PCone(cl);

}

/*
JUSTIFICATION
     LEFT     Align left edge of text to origin
     CENTER   Align center of text to origin
     RIGHT    Align right edge of text to origin

FILE FORMAT/DEFAULTS
     AsciiText {
          string         ""    # MFString
          spacing        1     # SFFloat
          justification  LEFT  # SFEnum
          width          0     # MFFloat
     }
*/

int ParseVRML::PAsciiText(CAsciiText& C) {
        int token;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
                case TOKEN_STRING:
                        C.string=t->GetStringBB('"');
                        if(!C.string) t->Error(syntax_error);
                        break;
                case TOKEN_SPACING:
                        C.spacing=t->GetDouble();
                        break;
                case TOKEN_JUSTIFICATION:
                        switch(GetToken()) {
                        case TOKEN_LEFT:        C.justification=JUSTIFICATION_LEFT;break;
                        case TOKEN_RIGHT:       C.justification=JUSTIFICATION_RIGHT;break;
                        case TOKEN_CENTER:      C.justification=JUSTIFICATION_CENTER;break;
                        default:                        t->Error("Expected LEFT|RIGHT|CENTER");
                        }
                        break;
                case TOKEN_WIDTH:
                        C.width=t->GetDouble();
                        break;
                case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
PARTS
     SIDES       The conical part
     BOTTOM      The bottom circular face
     ALL         All parts

FILE FORMAT/DEFAULTS
     Cone {
          parts         ALL     # SFBitMask
          bottomRadius  1       # SFFloat
          height        2       # SFFloat
     }
*/

int ParseVRML::PCone(CCone& C) {
        int token,token1;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
                case TOKEN_PARTS:
                        C.parts=PARTS_CONE_NONE;
                        do {
                                switch(token1=GetToken()) {
                                case TOKEN_SIDES:       C.parts|=PARTS_CONE_SIDES;break;
                                case TOKEN_BOTTOM:      C.parts|=PARTS_CONE_BOTTOM;break;
                                case TOKEN_ALL:         C.parts|=PARTS_CONE_ALL;break;
                                case TOKEN_PIPE:        break;
                                default:                        t->Error("Expected SIDES|BOTTOM|ALL");
                                } 
                        } while(token1==TOKEN_PIPE || ShowToken()==TOKEN_PIPE);
                        break;
                case TOKEN_BOTTOMRADIUS:
                        C.bottomRadius=t->GetDouble();
                        break;
                case TOKEN_HEIGHT:
                        C.height=t->GetDouble();
                        break;
                case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
FILE FORMAT/DEFAULTS
     Coordinate3 {
          point  0 0 0  # MFVec3f
     }
*/

int ParseVRML::PCoordinate3(CCoordinate3& C) {
        if(GetToken()!=TOKEN_LCB)   t->Error(syntax_error);
        if(GetToken()!=TOKEN_POINT) t->Error(syntax_error);
        PMFVec3f(C.size,&C.data);
        if(GetToken()!=TOKEN_RCB)   t->Error(syntax_error);
        return 0;
}

int ParseVRML::PMFVec3f(int& count,Vector3** data) {
        int token;
        token=ShowToken();
        if(token==TOKEN_NUMBER) {
                if(*data) delete *data;
                *data=new Vector3[1];
                (*data)[0].x=t->GetDouble();
                (*data)[0].y=t->GetDouble();
                (*data)[0].z=t->GetDouble();
                count=1;
                if(GetToken()!=TOKEN_RCB) t->Error(syntax_error);
                return 1;
        }
        if(token==TOKEN_LEB) {
                GetToken();
                NS_List<Vector3> list;
                NS_List<Vector3>* p;
                count=0;
                do {
                        p=new NS_List<Vector3>;
                        p->data.x=t->GetDouble();
                        p->data.y=t->GetDouble();
                        p->data.z=t->GetDouble();
                        count++;
                        list.InsertFirst(p);
                        token=GetToken();
                        if(token!=TOKEN_REB && token!=TOKEN_COMMA) 
                                t->Error("Cannot parse array");
                } while(token!=TOKEN_REB && ShowToken()!=TOKEN_REB);
                if(ShowToken()==TOKEN_REB) GetToken();
                p=list.GetNext();
                if(*data) delete *data;
                *data=new Vector3[count];
                int c=count;
                count=0;
                while(p) {
                        (*data)[c-count-1]=p->data;
                        p=p->GetNext();
                        count++;
                }
                return 1;       
        }
        t->Error(syntax_error);
        return 0;
}

int ParseVRML::PMFLong(int& count,long** data) {
        int token;
        token=ShowToken();
        if(token==TOKEN_NUMBER) {
                if(*data) delete *data;
                *data=new long[1];
                (*data)[0]=t->GetLong();
                count=1;
                if(GetToken()!=TOKEN_RCB) t->Error(syntax_error);
                return 1;
        }
        if(token==TOKEN_LEB) {
                GetToken();
                NS_List<long> list;
                NS_List<long>* p;
                count=0;
                do {
                        p=new NS_List<long>;
                        p->data=t->GetLong();
                        count++;
                        list.InsertFirst(p);
                        token=GetToken();
                        if(token!=TOKEN_REB && token!=TOKEN_COMMA) 
                                t->Error("Cannot parse array");
                } while(token!=TOKEN_REB && ShowToken()!=TOKEN_REB);
                if(ShowToken()==TOKEN_REB) GetToken();
                p=list.GetNext();
                if(*data) delete *data;
                *data=new long[count];
                int c=count;
                count=0;
                while(p) {
                        (*data)[c-count-1]=p->data;
                        p=p->GetNext();
                        count++;
                }
                return 1;       
        }
        t->Error(syntax_error);
        return 0;
}

int ParseVRML::PMFFloat(int& count,double** data) {
        int token;
        token=ShowToken();
        if(token==TOKEN_NUMBER) {
                if(*data) delete *data;
                *data=new double[1];
                (*data)[0]=t->GetDouble();
                count=1;
                if(GetToken()!=TOKEN_RCB) t->Error(syntax_error);
                return 1;
        }
        if(token==TOKEN_LEB) {
                GetToken();
                NS_List<double> list;
                NS_List<double>* p;
                count=0;
                do {
                        p=new NS_List<double>;
                        p->data=t->GetDouble();
                        count++;
                        list.InsertFirst(p);
                        token=GetToken();
                        if(token!=TOKEN_REB && token!=TOKEN_COMMA) 
                                t->Error("Cannot parse array");
                } while(token!=TOKEN_REB && ShowToken()!=TOKEN_REB);
                if(ShowToken()==TOKEN_REB) GetToken();
                p=list.GetNext();
                if(*data) delete *data;
                *data=new double[count];
                int c=count;
                count=0;
                while(p) {
                        (*data)[c-count-1]=p->data;
                        p=p->GetNext();
                        count++;
                }
                return 1;       
        }
        t->Error(syntax_error);
        return 0;
}

        
/*
FILE FORMAT/DEFAULTS
     Cube {
          width   2     # SFFloat
          height  2     # SFFloat
          depth   2     # SFFloat
     }
*/

int ParseVRML::PCube(CCube& C) {
        int token;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
                case TOKEN_DEPTH:
                        C.depth=t->GetDouble();
                        break;
                case TOKEN_WIDTH:
                        C.width=t->GetDouble();
                        break;
                case TOKEN_HEIGHT:
                        C.height=t->GetDouble();
                        break;
                case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
PARTS
     SIDES   The cylindrical part
     TOP     The top circular face
     BOTTOM  The bottom circular face
     ALL     All parts

FILE FORMAT/DEFAULTS
     Cylinder {
          parts   ALL   # SFBitMask
          radius  1     # SFFloat
          height  2     # SFFloat
     }
*/

int ParseVRML::PCylinder(CCylinder& C) {
        int token,token1;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {

                switch(token=GetToken()) {
                case TOKEN_PARTS:
                        C.parts=PARTS_CYLINDER_NONE;;
                        do {
                                switch(token1=GetToken()) {
                                case TOKEN_SIDES:       C.parts|=PARTS_CYLINDER_SIDES;break;
                                case TOKEN_BOTTOM:      C.parts|=PARTS_CYLINDER_BOTTOM;break;
                                case TOKEN_TOP:         C.parts|=PARTS_CYLINDER_TOP;break;
                                case TOKEN_ALL:         C.parts|=PARTS_CYLINDER_ALL;break;
                                case TOKEN_PIPE:        break;
                                default:                        t->Error("Expected SIDES|BOTTOM|TOP|ALL");
                                } 
                        } while(token1==TOKEN_PIPE || ShowToken()==TOKEN_PIPE);
                        break;
                case TOKEN_RADIUS:
                        C.radius=t->GetDouble();
                        break;
                case TOKEN_HEIGHT:
                        C.height=t->GetDouble();
                        break;
                case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
FILE FORMAT/DEFAULTS
     DirectionalLight {
          on         TRUE       # SFBool
          intensity  1          # SFFloat
          color      1 1 1      # SFColor
          direction  0 0 -1     # SFVec3f
     }
*/

int     ParseVRML::PDirectionalLight(CDirectionalLight& C) {
        int token;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
                case TOKEN_INTENSITY:
                        C.intensity=t->GetDouble();
                        break;
                case TOKEN_COLOR:
                        C.color[0]=t->GetDouble();
                        C.color[1]=t->GetDouble();
                        C.color[2]=t->GetDouble();
                        break;
                case TOKEN_DIRECTION:
                        C.direction.x=t->GetDouble();
                        C.direction.y=t->GetDouble();
                        C.direction.z=t->GetDouble();
                        break;
                case TOKEN_ON:
                        switch(GetToken()) {
                        case TOKEN_TRUE:        C.on=1;break;
                        case TOKEN_FALSE:       C.on=0;break;
                        default:                        t->Error(syntax_error);
                        }
                        break;
                case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
     SERIF       Serif style (such as TimesRoman)
     SANS        Sans Serif Style (such as Helvetica)
     TYPEWRITER  Fixed pitch style (such as Courier)

STYLE
     NONE        No modifications to family
     BOLD        Embolden family
     ITALIC      Italicize or Slant family

FILE FORMAT/DEFAULTS
     FontStyle {
          size     10      # SFFloat
          family   SERIF   # SFEnum
          style    NONE    # SFBitMask
     }
*/

int ParseVRML::PFontStyle(CFontStyle& C) {
        int token,token1;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
                case TOKEN_PARTS:
                        C.style=STYLE_NONE;
                        do {
                                switch(token1=GetToken()) {
                                case TOKEN_NONE:        C.style=STYLE_NONE;break;
                                case TOKEN_BOLD:        C.style|=STYLE_BOLD;break;
                                case TOKEN_ITALIC:      C.style|=STYLE_ITALIC;break;
        case TOKEN_PIPE:        break;
                                default:                        t->Error("Expected SIDES|BOTTOM|TOP|ALL");
                                } 
                        } while(token1==TOKEN_PIPE || ShowToken()==TOKEN_PIPE);
                        break;
                case TOKEN_FAMILY:
                        switch(GetToken()) {
                        case TOKEN_SERIF:       C.family=FAMILY_SERIF;break;
                        case TOKEN_SANS:        C.family=FAMILY_SANS;break;
                        case TOKEN_TYPEWRITER: C.family=FAMILY_TYPEWRITER;break;
        default:                        t->Error(syntax_error);
                        }
                        break;
                case TOKEN_SIZE:
                        C.size=t->GetDouble();
                        break;
                case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
FILE FORMAT/DEFAULTS
     IndexedFaceSet {
          coordIndex         0  # MFLong
          materialIndex      -1 # MFLong
          normalIndex        -1 # MFLong
          textureCoordIndex  -1 # MFLong
     }
*/

int ParseVRML::PIndexedFaceSet(CIndexedFaceSet& C) {
        int token;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
                case TOKEN_COORDINDEX:          PMFLong(C.coordIndex_size,&C.coordIndex);break;
                case TOKEN_MATERIALINDEX:       PMFLong(C.materialIndex_size,&C.materialIndex);break;
                case TOKEN_NORMALINDEX:         PMFLong(C.normalIndex_size,&C.normalIndex);break;
                case TOKEN_TEXTURECOORDINDEX: PMFLong(C.textureCoordIndex_size,&C.textureCoordIndex);break;
                case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
FILE FORMAT/DEFAULTS
     IndexedLineSet {
          coordIndex         0  # MFLong
          materialIndex      -1 # MFLong
          normalIndex        -1 # MFLong
          textureCoordIndex  -1 # MFLong
     }
*/

int ParseVRML::PIndexedLineSet(CIndexedLineSet& C) {
        int token;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
                case TOKEN_COORDINDEX:          PMFLong(C.coordIndex_size,&C.coordIndex);break;
                case TOKEN_MATERIALINDEX:       PMFLong(C.materialIndex_size,&C.materialIndex);break;
                case TOKEN_NORMALINDEX:         PMFLong(C.normalIndex_size,&C.normalIndex);break;
                case TOKEN_TEXTURECOORDINDEX: PMFLong(C.textureCoordIndex_size,&C.textureCoordIndex);break;
                case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
    Info {
          string  "<Undefined info>"      # SFString
     }
*/

int ParseVRML::PInfo(CInfo& C) {
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        if(GetToken()!=TOKEN_STRING) t->Error(syntax_error);
        if(C.string) delete C.string;
        C.string=t->GetStringBB('"');
        if(GetToken()!=TOKEN_RCB) t->Error(syntax_error);
        return 1;
}

/*
FILE FORMAT/DEFAULTS
     LOD {
          range [ ]    # MFFloat
          center 0 0 0  # SFVec3f
     }
*/

int ParseVRML::PLOD(CLOD& C) {
        int token;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
        case TOKEN_CENTERs:
            C.center[0]=t->GetDouble();
            C.center[1]=t->GetDouble();
            C.center[2]=t->GetDouble();
            break;
        case TOKEN_RANGE:
            return PMFFloat(C.range_count,&C.range);
        case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
FILE FORMAT/DEFAULTS
     Material {
          ambientColor   0.2 0.2 0.2    # MFColor
          diffuseColor   0.8 0.8 0.8    # MFColor
          specularColor  0 0 0          # MFColor
          emissiveColor  0 0 0          # MFColor
          shininess      0.2            # MFFloat
          transparency   0              # MFFloat
     }
*/

int ParseVRML::PMaterial(CMaterial& C) {
        int token;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
        case TOKEN_AMBIENTCOLOR:
            C.ambientColor[0]=t->GetDouble();
            C.ambientColor[1]=t->GetDouble();
            C.ambientColor[2]=t->GetDouble();
            break;
        case TOKEN_DIFFUSECOLOR:
            C.diffuseColor[0]=t->GetDouble();
            C.diffuseColor[1]=t->GetDouble();
            C.diffuseColor[2]=t->GetDouble();
            break;
        case TOKEN_SPECULARCOLOR:
            C.specularColor[0]=t->GetDouble();
            C.specularColor[1]=t->GetDouble();
            C.specularColor[2]=t->GetDouble();
            break;
        case TOKEN_EMISSIVECOLOR:
            C.emissiveColor[0]=t->GetDouble();
            C.emissiveColor[1]=t->GetDouble();
            C.emissiveColor[2]=t->GetDouble();
            break;
        case TOKEN_SHININESS:
            C.shininess=t->GetDouble();
            break;
        case TOKEN_TRANSPARENCY:
            C.transparency=t->GetDouble();
            break;
        case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
BINDINGS
     DEFAULT            Use default binding
     OVERALL            Whole object has same material
     PER_PART           One material for each part of object
     PER_PART_INDEXED   One material for each part, indexed
     PER_FACE           One material for each face of object
     PER_FACE_INDEXED   One material for each face, indexed
     PER_VERTEX         One material for each vertex of object
     PER_VERTEX_INDEXED One material for each vertex, indexed

FILE FORMAT/DEFAULTS
     MaterialBinding {
          value  DEFAULT        # SFEnum
     }
*/

int ParseVRML::PMaterialBinding(CMaterialBinding& C) {
        int token;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
        case TOKEN_VALUE: switch(GetToken()) {
            case TOKEN_DEFAULT: C.value=VALUE_DEFAULT;break;
            case TOKEN_OVERALL: C.value=VALUE_OVERALL;break;
            case TOKEN_PER_PART: C.value=VALUE_PER_PART;break;
            case TOKEN_PER_PART_INDEXED: C.value=VALUE_PER_PART_INDEXED;break;
            case TOKEN_PER_FACE: C.value=VALUE_PER_FACE;break;
            case TOKEN_PER_FACE_INDEXED: C.value=VALUE_PER_FACE_INDEXED;break;
            case TOKEN_PER_VERTEX: C.value=VALUE_PER_VERTEX;break;
            case TOKEN_PER_VERTEX_INDEXED: C.value=VALUE_PER_VERTEX_INDEXED;break;
            }
            break;
        case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
FILE FORMAT/DEFAULTS
     MatrixTransform {
          matrix  1 0 0 0       # SFMatrix
                  0 1 0 0
                  0 0 1 0
                  0 0 0 1
     }
*/

int ParseVRML::PMatrixTransform(CMatrixTransform& C) {
        int token;
    int i,j;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
        case TOKEN_MATRIX:
            for(i=0;i<4;i++) for(j=0;j<4;j++) C.matrix[i][j]=t->GetDouble();
            break;
        case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
FILE FORMAT/DEFAULTS
     Normal {
          vector  0 0 1 # MFVec3f
     }
*/

int ParseVRML::PNormal(CNormal& C) {
        int token;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
	    switch(token=GetToken()) {
	    case TOKEN_VECTOR:
		PMFVec3f(C.size,&C.data);
		/*
		  C.vector[0]=t->GetDouble();
		  C.vector[1]=t->GetDouble();
		  C.vector[2]=t->GetDouble();
		*/
		break;
	    case TOKEN_RCB:
		break;
	    default: 
		t->Error(syntax_error);
	    }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
BINDINGS
     DEFAULT            Use default binding
     OVERALL            Whole object has same normal
     PER_PART           One normal for each part of object
     PER_PART_INDEXED   One normal for each part, indexed
     PER_FACE           One normal for each face of object
     PER_FACE_INDEXED   One normal for each face, indexed
     PER_VERTEX         One normal for each vertex of object
     PER_VERTEX_INDEXED One normal for each vertex, indexed

FILE FORMAT/DEFAULTS
     NormalBinding {
          value  DEFAULT        # SFEnum
     }
*/

int ParseVRML::PNormalBinding(CNormalBinding& C) {
        int token;
        if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
        do {
                switch(token=GetToken()) {
        case TOKEN_VALUE: switch(GetToken()) {
            case TOKEN_DEFAULT: C.value=VALUE_DEFAULT;break;
            case TOKEN_OVERALL: C.value=VALUE_OVERALL;break;
            case TOKEN_PER_PART: C.value=VALUE_PER_PART;break;
            case TOKEN_PER_PART_INDEXED: C.value=VALUE_PER_PART_INDEXED;break;
            case TOKEN_PER_FACE: C.value=VALUE_PER_FACE;break;
            case TOKEN_PER_FACE_INDEXED: C.value=VALUE_PER_FACE_INDEXED;break;
            case TOKEN_PER_VERTEX: C.value=VALUE_PER_VERTEX;break;
            case TOKEN_PER_VERTEX_INDEXED: C.value=VALUE_PER_VERTEX_INDEXED;break;
            }
            break;
        case TOKEN_RCB:
                        break;
                default: 
                        t->Error(syntax_error);
                }
        } while(token!=TOKEN_RCB);
        return 0;
}

/*
FILE FORMAT/DEFAULTS
     OrthographicCamera {
          position         0 0 1        # SFVec3f
          orientation      0 0 1  0     # SFRotation
          focalDistance    5            # SFFloat
          height           2            # SFFloat
     }
*/

int ParseVRML::POrthographicCamera(COrthographicCamera& C) {
    int token;
    if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
    do {
        switch(token=GetToken()) {
        case TOKEN_POSITION:
            C.position[0]=t->GetDouble();
            C.position[1]=t->GetDouble();
            C.position[2]=t->GetDouble();
            break;
        case TOKEN_ORIENTATION:
            C.orientation[0]=t->GetDouble();
            C.orientation[1]=t->GetDouble();
            C.orientation[2]=t->GetDouble();
            C.orientation[3]=t->GetDouble();
            break;
        case TOKEN_FOCALDISTANCE:
            C.focalDistance=t->GetDouble();
            break;
        case TOKEN_HEIGHT:
            C.height=t->GetDouble();
            break;
        case TOKEN_RCB:
            break;
        default:
            break;
        }
    } while(token!=TOKEN_RCB);
    return 0;
}

/*
FILE FORMAT/DEFAULTS
     PerspectiveCamera {
          position         0 0 1        # SFVec3f
          orientation      0 0 1  0     # SFRotation
          focalDistance    5            # SFFloat
          heightAngle      0.785398     # SFFloat
     }
*/

int ParseVRML::PPerspectiveCamera(CPerspectiveCamera& C) {
    int token;
    if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
    do {
        switch(token=GetToken()) {
        case TOKEN_POSITION:
            C.position[0]=t->GetDouble();
            C.position[1]=t->GetDouble();
            C.position[2]=t->GetDouble();
            break;
        case TOKEN_ORIENTATION:
            C.orientation[0]=t->GetDouble();
            C.orientation[1]=t->GetDouble();
            C.orientation[2]=t->GetDouble();
            C.orientation[3]=t->GetDouble();
            break;
        case TOKEN_FOCALDISTANCE:
            C.focalDistance=t->GetDouble();
            break;
        case TOKEN_HEIGHTANGLE:
            C.heightAngle=t->GetDouble();
            break;
        case TOKEN_RCB:
            break;
        default: 
            t->Error(syntax_error);
        }
    } while(token!=TOKEN_RCB);
    return 0;
}

/*
FILE FORMAT/DEFAULTS
     PointLight {
          on         TRUE       # SFBool
          intensity  1          # SFFloat
          color      1 1 1      # SFColor
          location   0 0 1      # SFVec3f
     }
*/

int ParseVRML::PPointLight(CPointLight& C) {
    int token;
    if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
    do {
        switch(token=GetToken()) {
        case TOKEN_ON:
            switch(GetToken()) {
            case TOKEN_TRUE: C.on=1;break;
            case TOKEN_FALSE: C.on=0;break;
            default: t->Error(syntax_error);break;
            }
            break;
        case TOKEN_INTENSITY:
            C.intensity=t->GetDouble();
            break;
        case TOKEN_COLOR:
            C.color[0]=t->GetDouble();
            C.color[1]=t->GetDouble();
            C.color[2]=t->GetDouble();
            break;
        case TOKEN_LOCATION:
            C.location[0]=t->GetDouble();
            C.location[1]=t->GetDouble();
            C.location[2]=t->GetDouble();
            break;
        case TOKEN_RCB:
            break;
        default: 
            t->Error(syntax_error);
        }
    } while(token!=TOKEN_RCB);
    return 0;
}

/*
FILE FORMAT/DEFAULTS
     PointSet {
          startIndex  0 # SFLong
          numPoints   -1        # SFLong
     }
*/

int ParseVRML::PPointSet(CPointSet& C) {
    int token;
    if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
    do {
        switch(token=GetToken()) {
        case TOKEN_STARTINDEX:
            C.startIndex=t->GetLong();
            break;
        case TOKEN_NUMPOINTS:
            C.numPoints=t->GetLong();
            break;
        case TOKEN_RCB:
            break;
        default: 
            t->Error(syntax_error);
        }
    } while(token!=TOKEN_RCB);
    return 0;
}

/*
FILE FORMAT/DEFAULTS
     Rotation {
          rotation  0 0 1  0    # SFRotation
     }
*/

int ParseVRML::PRotation(CRotation& C) {
    int token;
    if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
    do {
        switch(token=GetToken()) {
        case TOKEN_ROTATIONs:
            C.axis.x=t->GetDouble();
            C.axis.y=t->GetDouble();
            C.axis.z=t->GetDouble();
            C.angle=t->GetDouble();
        case TOKEN_RCB:
            break;
        default: 
            t->Error(syntax_error);
        }
    } while(token!=TOKEN_RCB);
    return 0;
}

/*
FILE FORMAT/DEFAULTS
     Scale {
          scaleFactor  1 1 1    # SFVec3f
     }
*/

int ParseVRML::PScale(CScale& C) {
    int token;
    if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
    do {
        switch(token=GetToken()) {
        case TOKEN_SCALEFACTOR:
            C.scaleFactor[0]=t->GetDouble();
            C.scaleFactor[1]=t->GetDouble();
            C.scaleFactor[2]=t->GetDouble();
            break;
        case TOKEN_RCB:
            break;
        default: 
            t->Error(syntax_error);
        }
    } while(token!=TOKEN_RCB);
    return 0;
}

/*
VERTEX ORDERING ENUMS
     UNKNOWN_ORDERING    Ordering of vertices is unknown
     CLOCKWISE           Face vertices are ordered clockwise
                          (from the outside)
     COUNTERCLOCKWISE    Face vertices are ordered counterclockwise
                          (from the outside)

SHAPE TYPE ENUMS
     UNKNOWN_SHAPE_TYPE  Nothing is known about the shape
     SOLID               The shape encloses a volume

FACE TYPE ENUMS
     UNKNOWN_FACE_TYPE   Nothing is known about faces
     CONVEX              All faces are convex

FILE FORMAT/DEFAULTS
     ShapeHints {
          vertexOrdering  UNKNOWN_ORDERING      # SFEnum
          shapeType       UNKNOWN_SHAPE_TYPE    # SFEnum
          faceType        CONVEX                # SFEnum
          creaseAngle     0.5                   # SFFloat
     }
*/

int ParseVRML::PShapeHints(CShapeHints& C) {
    int token;
    if(GetToken()!=TOKEN_LCB) t->Error(syntax_error);
    do {
        switch(token=GetToken()) {
        case TOKEN_VERTEXORDERING:
            switch(GetToken()) {
            case TOKEN_UNKNOWN_ORDERING: C.vertexOrdering=VERTEXORDERING_UNKNOWN_ORDERING;break;
            case TOKEN_CLOCKWISE: C.vertexOrdering=VERTEXORDERING_CLOCKWISE;break;
            case TOKEN_COUNTERCLOCKWISE: C.vertexOrdering=VERTEXORDERING_COUNTERCLOCKWISE;break;
            default: t->Error(syntax_error);
            }
            break;
        case TOKEN_SHAPETYPE:
            switch(GetToken()) {
            case TOKEN_UNKNOWN_SHAPE_TYPE: C.shapeType=SHAPETYPE_UNKNOWN_SHAPE_TYPE;break;
            case TOKEN_SOLID: C.shapeType=SHAPETYPE_SOLID;break;
            default: t->Error(syntax_error);
            }
            break;
        case TOKEN_FACETYPE:
            switch(GetToken()) {
            case TOKEN_UNKNOWN_FACE_TYPE: C.faceType=FACETYPE_UNKNOWN_FACE_TYPE;break;
            case TOKEN_CONVEX: C.faceType=FACETYPE_CONVEX;break;
            default: t->Error(syntax_error);
            }
            break;
        case TOKEN_CREASEANGLE:
            C.creaseAngle=t->GetDouble();
            break;
        case TOKEN_RCB:
            break;
        default: 
            t->Error(syntax_error);
        }
    } while(token!=TOKEN_RCB);
    return 0;
}

// ***********************************************************************************

class Angle {
public:
    double angle,weight;
    inline bool operator>(const Angle& ang) {return (this->angle>ang.angle);}
    inline bool operator==(const Angle& ang) {return (this->angle==ang.angle);}
};



class Config {
public:
    int     rot;                    // optimize angle, bits 0-2 used
    double  angles[3];
    char*   outputvrml;
    char*   printpnt;
    char*   outputnav;
    double  bound[6];

            Config();
           ~Config();
    void    Parse(char* filename);
};

Config::Config() {
    int i;
    rot=0;
    outputvrml=0;
    printpnt=0;
    outputnav=0;
    angles[0]=angles[1]=angles[2]=0;
    for(i=0;i<6;i++) bound[i]=0;
}

Config::~Config() {
    if(outputvrml) delete outputvrml;
    if(printpnt) delete printpnt;
    if(outputnav) delete outputnav;
}

#define CTOKEN_COMMENT 4
#define CTOKEN_OUTPUTVRML 5
#define CTOKEN_OUTPUTPOINT 6
#define CTOKEN_OUTPUTNAV 7
#define CTOKEN_BOUND 11
#define CTOKEN_COMMA 12

Token cfg_token[]={
    {";",CTOKEN_COMMENT},
    {"outputvrml",CTOKEN_OUTPUTVRML},
    {"outputpnt",CTOKEN_OUTPUTPOINT},
    {"outputnav",CTOKEN_OUTPUTNAV},
    {"bound",CTOKEN_BOUND},
    {",",CTOKEN_COMMA},


    {0,0}};

void Config::Parse(char* filename) {
    Tokenizer t(filename);
    int token;
    int i;
    while((token=t.GetToken(cfg_token))!=TOKEN_EOF) {
        switch(token) {
        case CTOKEN_COMMENT: t.NewLine();break;
        case CTOKEN_OUTPUTVRML: outputvrml=t.GetString();break;
        case CTOKEN_OUTPUTPOINT: printpnt=t.GetString();break;
        case CTOKEN_OUTPUTNAV: outputnav=t.GetString();break;
        case CTOKEN_BOUND: bound[0]=t.GetDouble();
            for(i=1;i<6;i++) {
                if(t.GetToken(cfg_token)!=CTOKEN_COMMA) t.Error("comma expected.");
                bound[i]=t.GetDouble();
            }
            break;

        }
    }
}

void CIndexedFaceSet::CloseCoordIndex() {
    int i=0,j=0,k,count=0;
    while(j<coordIndex_size) {
        while(coordIndex[j]!=-1 && j<coordIndex_size) j++;
        if(coordIndex[i]!=coordIndex[j-1]) count++;
        i=j+1;j=i;
    }
    i=0;j=0;k=0;
    long* newcoordindex=new long[count+coordIndex_size];
    while(j<coordIndex_size) {
        while(coordIndex[j]!=-1) newcoordindex[k++]=coordIndex[j++];
        if(coordIndex[i]!=coordIndex[j-1]) newcoordindex[k++]=coordIndex[i];
        newcoordindex[k++]=-1;
        i=j+1;j=i;
    }
    delete coordIndex;
    coordIndex=newcoordindex;
    coordIndex_size+=count;
}


void ParseVRML::DoItQuickAndDirty() {
    int token;
    int i;
    int j,k,l,m,a,b,c;
    int open_brk=0,breakflag=0;

    int open_group=0, groupflag=0, groupcount=0, group_index=0;
    std::list<NS_List<CIndexedFaceSet>* >* SetListe;
    std::vector<std::list<NS_List<CIndexedFaceSet>* >* > GroupMap(1);

    
    FILE* f;

    NS_List<CCoordinate3>  clist;
    NS_List<CCoordinate3>* cp;
    NS_List<CIndexedFaceSet> flist;
    NS_List<CIndexedFaceSet>* fp;
    SortedList<double> xpos;
    SortedList<double> ypos;
    SortedList<double> zpos;
    SortedList<Angle> xang;
    SortedList<Angle> yang;
    SortedList<Angle> zang;

    t->SetEndTokenChars(etchars);
    if(t->GetToken(tok)!=TOKEN_HEAD) t->Error("Not a VRML 1.0 file.");
    t->NewLine();
    do {
        token=GetToken();
        switch(token) {
        case TOKEN_LCB:   
	    open_brk++;
	    breakflag=1;
	    if(groupflag) {
		open_group++;
		//fprintf(stdout, "group status open %d\n",open_group);
	    }
	    break;
        case TOKEN_RCB:     
	    open_brk--;
	    if(groupflag) {
		open_group--;
		//fprintf(stdout, "group status closed %d\n",open_group);
		if(open_group == 0) {
		    groupflag=0;
		    group_index++;
		    //fprintf(stdout, "Group closed!\n");
		}
	    }
	    break;
        case TOKEN_DEF:     delete t->GetString();break;
        case TOKEN_INFO:    {CInfo C;PInfo(C);} break;
	case TOKEN_ROTATION: {CRotation C;PRotation(C);} break;
	case TOKEN_GROUP: 
	    fprintf(stdout, "Detected one group.\n"); 
	    groupflag=1; groupcount++;
	    SetListe = new std::list<NS_List<CIndexedFaceSet>* >;
	    break;
        case TOKEN_COORDINATE3:  // relevant data for nav-file 
            cp=new NS_List<CCoordinate3>;
            PCoordinate3(cp->data);
            clist.InsertFirst(cp);
	    break;
        case TOKEN_INDEXEDFACESET: // relevant data for nav-file
            fp=new NS_List<CIndexedFaceSet>;
            PIndexedFaceSet(fp->data);
            flist.InsertFirst(fp);
	    if(groupflag) {
		GroupMap.push_back(SetListe);
		//fprintf(stdout, "This IndexSet belongs to a group. SetListe should exist...\n");
		SetListe->push_back(fp);
	    }
            break;
	case TOKEN_NORMAL: {CNormal C; PNormal(C);} break;
	case TOKEN_SHAPEHINTS: {CShapeHints C;PShapeHints(C);} break;
	case TOKEN_MATERIAL: {CMaterial C;PMaterial(C);} break;
	case TOKEN_NAME:    delete t->GetStringBB('"');break;
        }
    } while(open_brk>0 || breakflag==0);
    printf("parse ok.\n");
    

    
    Config cfg;
    cfg.Parse("vrml2nav.cfg");
    
    cp=clist.GetNext();                         // Szene nach (0,0,0)+bound setzen
    double xmin=1e20,ymin=1e20,zmin=1e20,xmax=-1e20,ymax=-1e20,zmax=-1e20;
    while(cp) {
        for(i=0;i<cp->data.size;i++) {
            if(xmin>cp->data.data[i].x) xmin=cp->data.data[i].x;
            if(ymin>cp->data.data[i].y) ymin=cp->data.data[i].y;
            if(zmin>cp->data.data[i].z) zmin=cp->data.data[i].z;
            if(xmax<cp->data.data[i].x) xmax=cp->data.data[i].x;
            if(ymax<cp->data.data[i].y) ymax=cp->data.data[i].y;
            if(zmax<cp->data.data[i].z) zmax=cp->data.data[i].z;
        }
        cp=cp->GetNext();
    }
    cp=clist.GetNext();
    while(cp) {
        for(i=0;i<cp->data.size;i++) {
            cp->data.data[i].x-=xmin-(xmax-xmin)*cfg.bound[0];
            cp->data.data[i].y-=ymin-(ymax-ymin)*cfg.bound[2];
            cp->data.data[i].z-=zmin-(zmax-zmin)*cfg.bound[4];
        }
                cp=cp->GetNext();
    }
    xmin=1e38;ymin=1e38;zmin=1e38;xmax=-1e38;ymax=-1e38;zmax=-1e38;
    cp=clist.GetNext();
    while(cp) {
        for(i=0;i<cp->data.size;i++) {
            if(xmin>cp->data.data[i].x) xmin=cp->data.data[i].x;
            if(ymin>cp->data.data[i].y) ymin=cp->data.data[i].y;
            if(zmin>cp->data.data[i].z) zmin=cp->data.data[i].z;
            if(xmax<cp->data.data[i].x) xmax=cp->data.data[i].x;
            if(ymax<cp->data.data[i].y) ymax=cp->data.data[i].y;
            if(zmax<cp->data.data[i].z) zmax=cp->data.data[i].z;
        }
        cp=cp->GetNext();
    }

    group_index=1; //FIXME: wg. push_back

    fp=flist.GetNext();
    if(cfg.outputnav) {
        f=fopen(cfg.outputnav,"w");
        if(!f) {printf("cannot create file %s\n",cfg.outputnav);exit(1);}
        fprintf(f,"dimension {\n  length <%f,%f,%f>\n  \
                   resolution <%i,%i,%i>\n}\n\n",
		xmax+(xmax-xmin)*(cfg.bound[0]+cfg.bound[1]),
		ymax+(ymax-ymin)*(cfg.bound[2]+cfg.bound[3]),
                zmax+(zmax-zmin)*(cfg.bound[4]+cfg.bound[5]),
		(int)((xmax/(xmax+ymax+zmax))*100),
		(int)((ymax/(xmax+ymax+zmax))*100),
		(int)((zmax/(xmax+ymax+zmax))*100));
        fprintf(f,"parameter {\n  prstep 50\n  itermax 100\n  Tfin 100.0\n  eps 1e-3\n  omega 1.8\n"
                  "  alpha 0.8\n  reynolds 100\n}\n");
        cp=clist.GetNext();
        while(cp && fp) {
	    if(groupcount==0) {
		fprintf(f,"poly {\n  points %d,\n",cp->data.size);
		for(i=0;i<cp->data.size;i++) {
		    fprintf(f,"    <%13.7f,%13.7f,%13.7f>%s\n",
			    cp->data.data[i].x, cp->data.data[i].y, cp->data.data[i].z,(i==cp->data.size-1)?"":",");
		}
		fp->data.CloseCoordIndex();
		fprintf(f,"vertices %d,\n    ",fp->data.coordIndex_size);
		for(i=0;i<fp->data.coordIndex_size;i++) 
		    fprintf(f,"%6d%s%s",fp->data.coordIndex[i],(i==fp->data.coordIndex_size-1)?"":",",
			    (fp->data.coordIndex[i]==-1)?"\n    ":"");
		fprintf(f,"}\n");
		cp=cp->GetNext();
		fp=fp->GetNext();
	    } else {
		fprintf(f,"poly {\n  points %d,\n",cp->data.size);
		for(i=0;i<cp->data.size;i++) {
		    fprintf(f,"    <%13.7f,%13.7f,%13.7f>%s\n",
			    cp->data.data[i].x, cp->data.data[i].y, cp->data.data[i].z,(i==cp->data.size-1)?"":",");
		}
		
		SetListe = GroupMap[group_index];
		int num_vertices=0;
		for (std::list<NS_List<CIndexedFaceSet>* >::iterator it=SetListe->begin();it!=SetListe->end();it++) {
		    (*it)->data.CloseCoordIndex();
		    num_vertices += (*it)->data.coordIndex_size;
		}
		fprintf(f,"vertices %d,\n    ",num_vertices);
		for (std::list<NS_List<CIndexedFaceSet>* >::iterator it=SetListe->begin();it!=SetListe->end();it++) {
		    (*it)->data.CloseCoordIndex();
		    for(i=0;i<(*it)->data.coordIndex_size;i++) {
			fprintf(f,"%6d ",(*it)->data.coordIndex[i]);
			if( (i==(*it)->data.coordIndex_size-1) && ( (*it)==SetListe->back() ))
			    fprintf(f, "\n");
			else
			    fprintf(f, ", ");
			if( (*it)->data.coordIndex[i]==-1) fprintf(f, "\n");


			    
		    }
		}
		fprintf(f,"}\n");
		cp=cp->GetNext();
		group_index++;
	    }
	}
	fclose(f);
	
    }
    group_index=1; //FIXME: wg. push_back
    if(cfg.outputvrml) {
        f=fopen(cfg.outputvrml,"w");
        if(!f) {printf("cannot create file %s\n",cfg.outputvrml);exit(1);}
        printf("Writing VRML file %s\n",cfg.outputvrml);
        fprintf(f,"#VRML V1.0 ascii\nSeparator{\n");
        cp=clist.GetNext(); 
        fp=flist.GetNext();
        while(cp) {
            fprintf(f,"Separator {\nCoordinate3{\npoint[");
            for(i=0;i<cp->data.size;i++) fprintf(f,"%f %f %f,",cp->data.data[i].x,cp->data.data[i].y,cp->data.data[i].z);
	    if(groupcount==0) {
		fprintf(f,"]\n}\nIndexedFaceSet {\ncoordIndex[");
		for(i=0;i<fp->data.coordIndex_size;i++) fprintf(f,"%d,",fp->data.coordIndex[i]);
		fprintf(f,"]\n}\n}\n");
		cp=cp->GetNext();
		fp=fp->GetNext();
	    } else {
		//FIXME:
		SetListe = GroupMap[group_index];
		fprintf(stdout, "Writing Group with %d sets\n", SetListe->size());
		fprintf(f,"]\n}\nGroup {\n");
		
		for (std::list<NS_List<CIndexedFaceSet>* >::iterator it=SetListe->begin();it!=SetListe->end();it++) {
		    //fprintf(stdout, "Writing Set\n");
		    fprintf(f, "Separator {\n");
		    fprintf(f, "IndexedFaceSet {\ncoordIndex[");
		    for(i=0;i<(*it)->data.coordIndex_size;i++) fprintf(f,"%d,",(*it)->data.coordIndex[i]);
		    fprintf(f, "]\n}\n}\n");
		}
		fprintf(f,"\n}\n}\n");
		cp=cp->GetNext();
		group_index++;
	    }
        }
        fprintf(f,"}");
        fclose(f);
    }
    if(cfg.printpnt) {
        f=fopen(cfg.printpnt,"w");
        if(!f) {printf("cannot create file %s\n",cfg.printpnt);exit(1);}
        printf("Writing point info %s\n",cfg.printpnt);
        SortedList<double>* p;
        p=xpos.GetNext();
        double last=0;
        while(p) {
            fprintf(f,"%15.10f %f\n",p->data,fabs(last-p->data));last=p->data;
            p=p->GetNext();
        }
        last=0;
        p=ypos.GetNext();
        fprintf(f,"\n");
        while(p) {
            fprintf(f,"%15.10f %f\n",p->data,fabs(last-p->data));last=p->data;
            p=p->GetNext();
        }
        last=0;
        p=zpos.GetNext();
        fprintf(f,"\n");
        while(p) {
            fprintf(f,"%15.10f %f\n",p->data,fabs(last-p->data));last=p->data;
            p=p->GetNext();
        }
        fclose(f);
    }
}

int main(int argc,char* argv[]) {
    ParseVRML p(argv[1]);
    p.DoItQuickAndDirty();
    return(1);
}
