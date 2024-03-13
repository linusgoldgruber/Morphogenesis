
// ** NOTES **
// Vector code CANNOT be inlined in header file because of dependencies
//    across vector classes (error generated: "Use of undeclared class..")
//

//#include <stdlib.h>

#ifndef VECTOR_DEF
	#define VECTOR_DEF

	//#include "common_defs.h" // pasted below

	#ifndef COMMON_DEF
        #define COMMON_DEF

        // Global defs

        #define		USE_SHADOWS			// Disable if you don't have FBOs, or to greatly speed up demo.

        // #define		USE_JPEG

        //#define		BUILD_CL			// CL - Visual Studio 2005 only (as of April 2009)

        #define TEX_SIZE		2048
        #define LIGHT_NEAR		0.5
        #define LIGHT_FAR		300.0
        #define DEGtoRAD		(3.141592/180.0)


        #ifdef _MSC_VER
            //#include <windows.h>
        #else
            typedef unsigned int		DWORD;
        #endif
        typedef unsigned int			uint;
        typedef unsigned short int		ushort;

        #define CLRVAL			uint
        #define COLOR(r,g,b)	( (uint(r*255.0f)<<24) | (uint(g*255.0f)<<16) | (uint(b*255.0f)<<8) )
        #define COLORA(r,g,b,a)	( (uint(a*255.0f)<<24) | (uint(b*255.0f)<<16) | (uint(g*255.0f)<<8) | uint(r*255.0f) )
        #define ALPH(c)			(float((c>>24) & 0xFF)/255.0)
        #define BLUE(c)			(float((c>>16) & 0xFF)/255.0)
        #define GRN(c)			(float((c>>8)  & 0xFF)/255.0)
        #define RED(c)			(float( c      & 0xFF)/255.0)

    #endif


	//#include <stdlib.h>
	#define VECTOR_INITIALIZE				// Initializes vectors

	//class Vector2DC;						// Forward Referencing
	typedef struct {
		unsigned char x, y;
	} Vector2DC;
	//class Vector2DI;
	typedef struct {
		int x, y;
	} Vector2DI;
	//class Vector2DF;
	typedef struct {
		float x, y;
	} Vector2DF;
	//class Vector3DC;
	typedef struct {
		unsigned char x, y, z;
	} Vector3DC;
	//class Vector3DI;
	typedef struct {
		int x, y, z;
	} Vector3DI;
	//class Vector3DF;
	typedef struct {
		float x, y, z;
	} Vector3DF;
	//class Vector4DC;
	typedef struct {
		unsigned char x, y, z;
	} Vector4DC;
	//class Vector4DF;
	typedef struct {
		float x, y, z, w;
	} Vector4DF;
	//class MatrixF;
	//class Matrix4F;
/*
	// Vector2DC Declaration

	#define VNAME		2DC
	#define VTYPE		unsigned char

	#define LUNA_CORE

	class LUNA_CORE Vector2DC {
	public:
		VTYPE x, y;

		// Constructors/Destructors
		inline Vector2DC();
		inline ~Vector2DC();
		inline Vector2DC (VTYPE xa, VTYPE ya);
		inline Vector2DC (Vector2DC &op);
		inline Vector2DC (Vector2DI &op);
		inline Vector2DC (Vector2DF &op);
		inline Vector2DC (Vector3DC &op);
		inline Vector2DC (Vector3DI &op);
		inline Vector2DC (Vector3DF &op);
		inline Vector2DC (Vector4DF &op);

		inline Vector2DC &Set (VTYPE xa, VTYPE ya)	{x=xa; y=ya; return *this;}

		// Member Functions
		inline Vector2DC &operator= (Vector2DC &op);
		inline Vector2DC &operator= (Vector2DI &op);
		inline Vector2DC &operator= (Vector2DF &op);
		inline Vector2DC &operator= (Vector3DC &op);
		inline Vector2DC &operator= (Vector3DI &op);
		inline Vector2DC &operator= (Vector3DF &op);
		inline Vector2DC &operator= (Vector4DF &op);

		inline Vector2DC &operator+= (Vector2DC &op);
		inline Vector2DC &operator+= (Vector2DI &op);
		inline Vector2DC &operator+= (Vector2DF &op);
		inline Vector2DC &operator+= (Vector3DC &op);
		inline Vector2DC &operator+= (Vector3DI &op);
		inline Vector2DC &operator+= (Vector3DF &op);
		inline Vector2DC &operator+= (Vector4DF &op);

		inline Vector2DC &operator-= (Vector2DC &op);
		inline Vector2DC &operator-= (Vector2DI &op);
		inline Vector2DC &operator-= (Vector2DF &op);
		inline Vector2DC &operator-= (Vector3DC &op);
		inline Vector2DC &operator-= (Vector3DI &op);
		inline Vector2DC &operator-= (Vector3DF &op);
		inline Vector2DC &operator-= (Vector4DF &op);

		inline Vector2DC &operator*= (Vector2DC &op);
		inline Vector2DC &operator*= (Vector2DI &op);
		inline Vector2DC &operator*= (Vector2DF &op);
		inline Vector2DC &operator*= (Vector3DC &op);
		inline Vector2DC &operator*= (Vector3DI &op);
		inline Vector2DC &operator*= (Vector3DF &op);
		inline Vector2DC &operator*= (Vector4DF &op);

		inline Vector2DC &operator/= (Vector2DC &op);
		inline Vector2DC &operator/= (Vector2DI &op);
		inline Vector2DC &operator/= (Vector2DF &op);
		inline Vector2DC &operator/= (Vector3DC &op);
		inline Vector2DC &operator/= (Vector3DI &op);
		inline Vector2DC &operator/= (Vector3DF &op);
		inline Vector2DC &operator/= (Vector4DF &op);

		// Note: Cross product does not exist for 2D vectors (only 3D)

		inline double Dot(Vector2DC &v);
		inline double Dot(Vector2DI &v);
		inline double Dot(Vector2DF &v);

		inline double Dist (Vector2DC &v);
		inline double Dist (Vector2DI &v);
		inline double Dist (Vector2DF &v);
		inline double Dist (Vector3DC &v);
		inline double Dist (Vector3DI &v);
		inline double Dist (Vector3DF &v);
		inline double Dist (Vector4DF &v);

		inline double DistSq (Vector2DC &v);
		inline double DistSq (Vector2DI &v);
		inline double DistSq (Vector2DF &v);
		inline double DistSq (Vector3DC &v);
		inline double DistSq (Vector3DI &v);
		inline double DistSq (Vector3DF &v);
		inline double DistSq (Vector4DF &v);

		inline Vector2DC &Normalize (void);
		inline double Length (void);
		inline VTYPE *Data (void);
	};

	#undef VNAME
	#undef VTYPE

	// Vector2DI Declaration

	#define VNAME		2DI
	#define VTYPE		int

	class LUNA_CORE Vector2DI {
	public:
		VTYPE x, y;

		// Constructors/Destructors
		inline Vector2DI();
		inline ~Vector2DI();
		inline Vector2DI (const VTYPE xa, const VTYPE ya);
		inline Vector2DI (const Vector2DC &op);
		inline Vector2DI (const Vector2DI &op);				// *** THESE SHOULD ALL BE const
		inline Vector2DI (const Vector2DF &op);
		inline Vector2DI (const Vector3DC &op);
		inline Vector2DI (const Vector3DI &op);
		inline Vector2DI (const Vector3DF &op);
		inline Vector2DI (const Vector4DF &op);

		// Member Functions
		inline Vector2DI &operator= (const Vector2DC &op);
		inline Vector2DI &operator= (const Vector2DI &op);
		inline Vector2DI &operator= (const Vector2DF &op);
		inline Vector2DI &operator= (const Vector3DC &op);
		inline Vector2DI &operator= (const Vector3DI &op);
		inline Vector2DI &operator= (const Vector3DF &op);
		inline Vector2DI &operator= (const Vector4DF &op);

		inline Vector2DI &operator+= (const Vector2DC &op);
		inline Vector2DI &operator+= (const Vector2DI &op);
		inline Vector2DI &operator+= (const Vector2DF &op);
		inline Vector2DI &operator+= (const Vector3DC &op);
		inline Vector2DI &operator+= (const Vector3DI &op);
		inline Vector2DI &operator+= (const Vector3DF &op);
		inline Vector2DI &operator+= (const Vector4DF &op);

		inline Vector2DI &operator-= (const Vector2DC &op);
		inline Vector2DI &operator-= (const Vector2DI &op);
		inline Vector2DI &operator-= (const Vector2DF &op);
		inline Vector2DI &operator-= (const Vector3DC &op);
		inline Vector2DI &operator-= (const Vector3DI &op);
		inline Vector2DI &operator-= (const Vector3DF &op);
		inline Vector2DI &operator-= (const Vector4DF &op);

		inline Vector2DI &operator*= (const Vector2DC &op);
		inline Vector2DI &operator*= (const Vector2DI &op);
		inline Vector2DI &operator*= (const Vector2DF &op);
		inline Vector2DI &operator*= (const Vector3DC &op);
		inline Vector2DI &operator*= (const Vector3DI &op);
		inline Vector2DI &operator*= (const Vector3DF &op);
		inline Vector2DI &operator*= (const Vector4DF &op);

		inline Vector2DI &operator/= (const Vector2DC &op);
		inline Vector2DI &operator/= (const Vector2DI &op);
		inline Vector2DI &operator/= (const Vector2DF &op);
		inline Vector2DI &operator/= (const Vector3DC &op);
		inline Vector2DI &operator/= (const Vector3DI &op);
		inline Vector2DI &operator/= (const Vector3DF &op);
		inline Vector2DI &operator/= (const Vector4DF &op);


		// Note: Cross product does not exist for 2D vectors (only 3D)

		inline double Dot (const Vector2DC &v);
		inline double Dot (const Vector2DI &v);
		inline double Dot (const Vector2DF &v);

		inline double Dist (const Vector2DC &v);
		inline double Dist (const Vector2DI &v);
		inline double Dist (const Vector2DF &v);
		inline double Dist (const Vector3DC &v);
		inline double Dist (const Vector3DI &v);
		inline double Dist (const Vector3DF &v);
		inline double Dist (const Vector4DF &v);

		inline double DistSq (const Vector2DC &v);
		inline double DistSq (const Vector2DI &v);
		inline double DistSq (const Vector2DF &v);
		inline double DistSq (const Vector3DC &v);
		inline double DistSq (const Vector3DI &v);
		inline double DistSq (const Vector3DF &v);
		inline double DistSq (const Vector4DF &v);

		inline Vector2DI &Normalize (void);
		inline double Length (void);
		inline VTYPE *Data (void);
	};

	#undef VNAME
	#undef VTYPE

	// Vector2DF Declarations

	#define VNAME		2DF
	#define VTYPE		float

	class LUNA_CORE Vector2DF {
	public:
		VTYPE x, y;

		// Constructors/Destructors
		 Vector2DF ();
		 ~Vector2DF ();
		 Vector2DF (const VTYPE xa, const VTYPE ya);
		 Vector2DF (const Vector2DC &op);
		 Vector2DF (const Vector2DI &op);
		 Vector2DF (const Vector2DF &op);
		 Vector2DF (const Vector3DC &op);
		 Vector2DF (const Vector3DI &op);
		 Vector2DF (const Vector3DF &op);
		 Vector2DF (const Vector4DF &op);

		 inline Vector2DF &Set (const float xa, const float ya);

		 // Member Functions
		 Vector2DF &operator= (const Vector2DC &op);
		 Vector2DF &operator= (const Vector2DI &op);
		 Vector2DF &operator= (const Vector2DF &op);
		 Vector2DF &operator= (const Vector3DC &op);
		 Vector2DF &operator= (const Vector3DI &op);
		 Vector2DF &operator= (const Vector3DF &op);
		 Vector2DF &operator= (const Vector4DF &op);

		 Vector2DF &operator+= (const Vector2DC &op);
		 Vector2DF &operator+= (const Vector2DI &op);
		 Vector2DF &operator+= (const Vector2DF &op);
		 Vector2DF &operator+= (const Vector3DC &op);
		 Vector2DF &operator+= (const Vector3DI &op);
		 Vector2DF &operator+= (const Vector3DF &op);
		 Vector2DF &operator+= (const Vector4DF &op);

		 Vector2DF &operator-= (const Vector2DC &op);
		 Vector2DF &operator-= (const Vector2DI &op);
		 Vector2DF &operator-= (const Vector2DF &op);
		 Vector2DF &operator-= (const Vector3DC &op);
		 Vector2DF &operator-= (const Vector3DI &op);
		 Vector2DF &operator-= (const Vector3DF &op);
		 Vector2DF &operator-= (const Vector4DF &op);

		 Vector2DF &operator*= (const Vector2DC &op);
		 Vector2DF &operator*= (const Vector2DI &op);
		 Vector2DF &operator*= (const Vector2DF &op);
		 Vector2DF &operator*= (const Vector3DC &op);
		 Vector2DF &operator*= (const Vector3DI &op);
		 Vector2DF &operator*= (const Vector3DF &op);
		 Vector2DF &operator*= (const Vector4DF &op);

		 Vector2DF &operator/= (const Vector2DC &op);
		 Vector2DF &operator/= (const Vector2DI &op);
		 Vector2DF &operator/= (const Vector2DF &op);
		 Vector2DF &operator/= (const Vector3DC &op);
		 Vector2DF &operator/= (const Vector3DI &op);
		 Vector2DF &operator/= (const Vector3DF &op);
		 Vector2DF &operator/= (const Vector4DF &op);

		 Vector2DF &operator/= (const double v)		{x /= (float) v; y /= (float) v; return *this;}

		// Note: Cross product does not exist for 2D vectors (only 3D)

		 double Dot(const Vector2DC &v);
		 double Dot(const Vector2DI &v);
		 double Dot(const Vector2DF &v);

		 double Dist (const Vector2DC &v);
		 double Dist (const Vector2DI &v);
		 double Dist (const Vector2DF &v);
		 double Dist (const Vector3DC &v);
		 double Dist (const Vector3DI &v);
		 double Dist (const Vector3DF &v);
		 double Dist (const Vector4DF &v);

		 double DistSq (const Vector2DC &v);
		 double DistSq (const Vector2DI &v);
		 double DistSq (const Vector2DF &v);
		 double DistSq (const Vector3DC &v);
		 double DistSq (const Vector3DI &v);
		 double DistSq (const Vector3DF &v);
		 double DistSq (const Vector4DF &v);

		 Vector2DF &Normalize (void);
		 double Length (void);
		 VTYPE *Data (void);
	};

	#undef VNAME
	#undef VTYPE

	// Vector3DC Declaration

	#define VNAME		3DC
	#define VTYPE		unsigned char

	class LUNA_CORE Vector3DC {
	public:
		VTYPE x, y, z;

		// Constructors/Destructors
		inline Vector3DC();
		inline ~Vector3DC();
		inline Vector3DC (const VTYPE xa, const VTYPE ya, const VTYPE za);
		inline Vector3DC  ( const Vector2DC &op);
		inline Vector3DC  ( const Vector2DI &op);
		inline Vector3DC  ( const Vector2DF &op);
		inline Vector3DC  ( const Vector3DC &op);
		inline Vector3DC  ( const Vector3DI &op);
		inline Vector3DC  ( const Vector3DF &op);
		inline Vector3DC  ( const Vector4DF &op);

		// Member Functions
		inline Vector3DC &Set (VTYPE xa, VTYPE ya, VTYPE za);

		inline Vector3DC &operator=  ( const Vector2DC &op);
		inline Vector3DC &operator=  ( const Vector2DI &op);
		inline Vector3DC &operator=  ( const Vector2DF &op);
		inline Vector3DC &operator=  ( const Vector3DC &op);
		inline Vector3DC &operator=  ( const Vector3DI &op);
		inline Vector3DC &operator=  ( const Vector3DF &op);
		inline Vector3DC &operator=  ( const Vector4DF &op);

		inline Vector3DC &operator+=  ( const Vector2DC &op);
		inline Vector3DC &operator+=  ( const Vector2DI &op);
		inline Vector3DC &operator+=  ( const Vector2DF &op);
		inline Vector3DC &operator+=  ( const Vector3DC &op);
		inline Vector3DC &operator+=  ( const Vector3DI &op);
		inline Vector3DC &operator+=  ( const Vector3DF &op);
		inline Vector3DC &operator+=  ( const Vector4DF &op);

		inline Vector3DC &operator-=  ( const Vector2DC &op);
		inline Vector3DC &operator-=  ( const Vector2DI &op);
		inline Vector3DC &operator-=  ( const Vector2DF &op);
		inline Vector3DC &operator-=  ( const Vector3DC &op);
		inline Vector3DC &operator-=  ( const Vector3DI &op);
		inline Vector3DC &operator-=  ( const Vector3DF &op);
		inline Vector3DC &operator-=  ( const Vector4DF &op);

		inline Vector3DC &operator*=  ( const Vector2DC &op);
		inline Vector3DC &operator*=  ( const Vector2DI &op);
		inline Vector3DC &operator*=  ( const Vector2DF &op);
		inline Vector3DC &operator*=  ( const Vector3DC &op);
		inline Vector3DC &operator*=  ( const Vector3DI &op);
		inline Vector3DC &operator*=  ( const Vector3DF &op);
		inline Vector3DC &operator*=  ( const Vector4DF &op);

		inline Vector3DC &operator/=  ( const Vector2DC &op);
		inline Vector3DC &operator/=  ( const Vector2DI &op);
		inline Vector3DC &operator/=  ( const Vector2DF &op);
		inline Vector3DC &operator/=  ( const Vector3DC &op);
		inline Vector3DC &operator/=  ( const Vector3DI &op);
		inline Vector3DC &operator/=  ( const Vector3DF &op);
		inline Vector3DC &operator/=  ( const Vector4DF &op);

		inline Vector3DC &Cross  ( const Vector3DC &v);
		inline Vector3DC &Cross  ( const Vector3DI &v);
		inline Vector3DC &Cross  ( const Vector3DF &v);

		inline double Dot ( const Vector3DC &v);
		inline double Dot ( const Vector3DI &v);
		inline double Dot ( const Vector3DF &v);

		inline double Dist  ( const Vector2DC &v);
		inline double Dist  ( const Vector2DI &v);
		inline double Dist  ( const Vector2DF &v);
		inline double Dist  ( const Vector3DC &v);
		inline double Dist  ( const Vector3DI &v);
		inline double Dist  ( const Vector3DF &v);
		inline double Dist  ( const Vector4DF &v);

		inline double DistSq  ( const Vector2DC &v);
		inline double DistSq  ( const Vector2DI &v);
		inline double DistSq  ( const Vector2DF &v);
		inline double DistSq  ( const Vector3DC &v);
		inline double DistSq  ( const Vector3DI &v);
		inline double DistSq  ( const Vector3DF &v);
		inline double DistSq  ( const Vector4DF &v);

		inline Vector3DC &Normalize (void);
		inline double Length (void);
		inline VTYPE *Data (void);
	};

	#undef VNAME
	#undef VTYPE
*/
	///////////Vector 3DF functions

    Vector3DI Vector3DI_init() {
        Vector3DI v;
        v.x = 0;
        v.y = 0;
        v.z = 0;
        return v;
    }

    void Vector3DI_destroy(Vector3DI *v) {
        // do nothing
    }

    Vector3DI Vector3DI_init_with_values(double xa, double ya, double za) {
        Vector3DI v;
        v.x = xa;
        v.y = ya;
        v.z = za;
        return v;
    }

    Vector3DI Vector3DI_init_with_Vector2DC(Vector2DI *op) {
        Vector3DI v;
        v.x = op->x;
        v.y = op->y;
        v.z = 0;
        return v;
    }

    Vector3DI Vector3DI_init_with_Vector2DI(Vector2DI *op) {
        Vector3DI v;
        v.x = op->x;
        v.y = op->y;
        v.z = 0;
        return v;
    }

    Vector3DI Vector3DI_init_with_Vector2DF(Vector2DF *op) {
        Vector3DI v;
        v.x = op->x;
        v.y = op->y;
        v.z = 0;
        return v;
    }

    Vector3DI Vector3DI_init_with_Vector3DC(Vector3DC *op) {
        Vector3DI v;
        v.x = op->x;
        v.y = op->y;
        v.z = op->z;
        return v;
    }

    Vector3DI Vector3DI_init_with_Vector3DI(Vector3DI *op) {
        Vector3DI v;
        v.x = op->x;
        v.y = op->y;
        v.z = op->z;
        return v;
    }

    Vector3DI Vector3DI_init_with_Vector3DF(Vector3DF *op) {
        Vector3DI v;
        v.x = op->x;
        v.y = op->y;
        v.z = op->z;
        return v;
    }

    Vector3DI Vector3DI_init_with_Vector4DF(Vector4DF *op) {
        Vector3DI v;
        v.x = op->x;
        v.y = op->y;
        v.z = op->z;
        return v;
    }

    Vector3DI *Vector3DI_Set(Vector3DI *v, double xa, double ya, double za) {
        v->x = xa;
        v->y = ya;
        v->z = za;
        return v;
    }

    Vector3DI *Vector3DI_operator_equal_int(Vector3DI *v, int op) {
        v->x = op;
        v->y = op;
        v->z = op;
        return v;
    }

    Vector3DI *Vector3DI_operator_equal_double(Vector3DI *v, double op) {
        v->x = op;
        v->y = op;
        v->z = op;
        return v;
    }

    Vector3DI *Vector3DI_operator_equal_Vector2DC(Vector3DI *v, Vector2DI *op) {
        v->x = op->x;
        v->y = op->y;
        v->z = 0;
        return v;
    }

    Vector3DI *Vector3DI_operator_equal_Vector2DI(Vector3DI *v, Vector2DI *op) {
        v->x = op->x;
        v->y = op->y;
        v->z = 0;
        return v;
    }

    Vector3DI *Vector3DI_operator_equal_Vector2DF(Vector3DI *v, Vector2DF *op) {
        v->x = op->x;
        v->y = op->y;
        v->z = 0;
        return v;
    }

    Vector3DI *Vector3DI_operator_equal_Vector3DC(Vector3DI *v, Vector3DC *op) {
        v->x = op->x;
        v->y = op->y;
        v->z = op->z;
        return v;
    }

    Vector3DI *Vector3DI_operator_equal_Vector3DI(Vector3DI *v, Vector3DI *op) {
        v->x = op->x;
        v->y = op->y;
        v->z = op->z;
        return v;
    }

//     Vector3DI *Vector3DI_operator_equal_Vector3DF(Vector3DF *v, Vector3DI *op) {
//         v->x = op->x;
//         v->y = op->y;
//         v->z = op->z;
//         return v;
//     }

    Vector3DI *Vector3DI_operator_equal_Vector4DF(Vector3DI *v, Vector4DF *op) {
        v->x = op->x;
        v->y = op->y;
        v->z = op->z;
        return v;
    }

    Vector3DF Vector3DI_to_Vector3DF(Vector3DI op) {
		Vector3DF result;
		result.x = (float) op.x;
		result.y = (float) op.y;
		result.z = (float) op.z;
		return result;
	}

    Vector3DI *Vector3DI_addInt(Vector3DI *vector, int op) {
		vector->x += op;
		vector->y += op;
		vector->z += op;
		return vector;
    }

    Vector3DI *Vector3DI_addDouble(Vector3DI *vector, double op) {
        vector->x += op;
        vector->y += op;
        vector->z += op;
        return vector;
    }

    Vector3DI *Vector3DI_addVector2DC(Vector3DI *vector, Vector2DI *op) {
        vector->x += op->x;
        vector->y += op->y;
        return vector;
    }

    Vector3DI *Vector3DI_addVector2DI(Vector3DI *vector, Vector2DI *op) {
        vector->x += op->x;
        vector->y += op->y;
        return vector;
    }

    Vector3DI *Vector3DI_addVector2DF(Vector3DI *vector, Vector2DF *op) {
        vector->x += op->x;
        vector->y += op->y;
        return vector;
    }

    Vector3DI *Vector3DI_addVector3DC(Vector3DI *vector, Vector3DC *op) {
        vector->x += op->x;
        vector->y += op->y;
        vector->z += op->z;
        return vector;
    }

    Vector3DI *Vector3DI_addVector3DI(Vector3DI *vector, Vector3DI *op) {
        vector->x += op->x;
        vector->y += op->y;
        vector->z += op->z;
        return vector;
    }

    Vector3DI *Vector3DI_addVector3DF(Vector3DI *vector, Vector3DF *op) {
        vector->x += (float) op->x;
        vector->y += (float) op->y;
        vector->z += (float) op->z;
        return vector;
    }

    Vector3DI *Vector3DI_addVector4DF(Vector3DI *vector, Vector4DF *op) {
        vector->x += op->x;
        vector->y += op->y;
        vector->z += op->z;
        return vector;
    }

    Vector3DI Vector3DI_subtractInt(Vector3DI *v, int op) {
        v->x -= op;
        v->y -= op;
        v->z -= op;
        return *v;
    }

    Vector3DI Vector3DI_subtractDouble(Vector3DI *v, double op) {
        v->x -= op;
        v->y -= op;
        v->z -= op;
        return *v;
    }

    Vector3DI Vector3DI_subtractVector2DC(Vector3DI *v, Vector2DI *op) {
        v->x -= op->x;
        v->y -= op->y;
        return *v;
    }

    Vector3DI Vector3DI_subtractVector2DI(Vector3DI *v, Vector2DI *op) {
        v->x -= op->x;
        v->y -= op->y;
        return *v;
    }

    Vector3DI Vector3DI_subtractVector2DF(Vector3DI *v, Vector2DF *op) {
        v->x -= op->x;
        v->y -= op->y;
        return *v;
    }

    Vector3DI Vector3DI_subtractVector3DC(Vector3DI *v, Vector3DC *op) {
        v->x -= op->x;
        v->y -= op->y;
        v->z -= op->z;
        return *v;
    }

    Vector3DI Vector3DI_subtractVector3DI(Vector3DI *v, Vector3DI *op) {
        v->x -= op->x;
        v->y -= op->y;
        v->z -= op->z;
        return *v;
    }

    Vector3DI Vector3DI_subtractVector3DF(Vector3DI *v, Vector3DF *op) {
        v->x -= (float) op->x;
        v->y -= (float) op->y;
        v->z -= (float) op->z;
        return *v;
    }

    Vector3DI Vector3DI_subtractVector4DF(Vector3DI *v, Vector4DF *op) {
        v->x -= op->x;
        v->y -= op->y;
        v->z -= op->z;
        return *v;
    }

    Vector3DI *Vector3DI_multiplyInt(Vector3DI *v, int *op) {
        v->x *= *op;
        v->y *= *op;
        v->z *= *op;
        return v;
    }

    Vector3DI *Vector3DI_multiplyDouble(Vector3DI *v, double op) {
        v->x *= op;
        v->y *= op;
        v->z *= op;
        return v;
    }

    Vector3DI *Vector3DI_multiplyVector2DC(Vector3DI *v, Vector2DI *op) {
        v->x *= op->x;
        v->y *= op->y;
    return v;
    }

    Vector3DI *Vector3DI_multiplyVector2DI(Vector3DI *v, Vector2DI *op) {
    v->x *= op->x;
    v->y *= op->y;
    return v;
}

    Vector3DI *Vector3DI_multiplyVector2DF(Vector3DI *v, Vector2DF *op) {
        v->x *= op->x;
        v->y *= op->y;
        return v;
    }

    Vector3DI *Vector3DI_multiplyVector3DC(Vector3DI *v, Vector3DC *op) {
        v->x *= op->x;
        v->y *= op->y;
        v->z *= op->z;
        return v;
    }

    Vector3DI *Vector3DI_multiplyVector3DI(Vector3DI *v, Vector3DI *op) {
        v->x *= op->x;
        v->y *= op->y;
        v->z *= op->z;
        return v;
    }
    Vector3DI *Vector3DI_multiplyVector3DF(Vector3DI *v, Vector3DF *op) {
        v->x *= (float) op->x;
        v->y *= (float) op->y;
        v->z *= (float) op->z;
        return v;
    }

    Vector3DI *Vector3DI_multiplyVector4DF(Vector3DI *v, Vector4DF *op) {
        v->x *= op->x;
        v->y *= op->y;
        v->z *= op->z;
        return v;
    }
/*

    Vector3DI *Vector3DI_multiplyMatrix4F(Vector3DI *v, Matrix4F *op) {
        // implementation MISSING
        return v;
    }

    Vector3DI *Vector3DI_multiplyMatrixF(Vector3DI *v, MatrixF *op) {
		// implementation MISSING
        return v;
    }
*/

    Vector3DI Vector3DI_divideByInt(Vector3DI *v, int op) {
        v->x /= op;
        v->y /= op;
        v->z /= op;
        return *v;
    }

    Vector3DI Vector3DI_divideByDouble(Vector3DI *v, double op) {
        v->x /= op;
        v->y /= op;
        v->z /= op;
        return *v;
    }

    Vector3DI Vector3DI_divideByVector2DC(Vector3DI *v, Vector2DC op) {
        v->x /= op.x;
        v->y /= op.y;
        return *v;
    }

    Vector3DI Vector3DI_divideByVector2DI(Vector3DI *v, Vector2DI op) {
        v->x /= op.x;
        v->y /= op.y;
        return *v;
    }

    Vector3DI Vector3DI_divideByVector2DF(Vector3DI *v, Vector2DF op) {
        v->x /= op.x;
        v->y /= op.y;
        return *v;
    }

    Vector3DI Vector3DI_divideByVector3DC(Vector3DI *v, Vector3DC op) {
        v->x /= op.x;
        v->y /= op.y;
        v->z /= op.z;
        return *v;
    }

    Vector3DI Vector3DI_divideByVector3DI(Vector3DI *v, Vector3DI op) {
        v->x /= op.x;
        v->y /= op.y;
        v->z /= op.z;
        return *v;
    }

    Vector3DI Vector3DI_divideByVector3DF(Vector3DI *v, Vector3DF op) {
        v->x /= (float) op.x;
        v->y /= (float) op.y;
        v->z /= (float) op.z;
        return *v;
    }

    Vector3DI Vector3DI_divideByVector4DF(Vector3DI *v, Vector4DF op) {
        v->x /= op.x;
        v->y /= op.y;
        v->z /= op.z;
        return *v;
    }
/*

    Vector3DI Vector3DI_addInt(Vector3DI *vector, int op) {

        Vector3DI result;

        result.x = vector->x + (float)op;

        result.y = vector->y + (float)op;

        result.z = vector->z + (float)op;

        return result;

    }*/


    Vector3DI Vector3DI_addFloat(Vector3DI *vector, float op) {
        Vector3DI result;
        result.x = vector->x + op;
        result.y = vector->y + op;
        result.z = vector->z + op;
        return result;
    }

    Vector3DI Vector3DI_addVector(Vector3DI *vector, Vector3DI *op) {
        Vector3DI result;
        result.x = vector->x + op->x;
        result.y = vector->y + op->y;
        result.z = vector->z + op->z;
        return result;
    }

    //Redefinition of Vector_3DF_subtractInt???

    /*Vector3DI Vector3DI_subtractInt(Vector3DI *vector, int op) {
        Vector3DI result;
        result.x = vector->x - (float)op;
        result.y = vector->y - (float)op;
        result.z = vector->z - (float)op;
        return result;
    }*/

    Vector3DI Vector3DI_subtractFloat(Vector3DI *vector, float op) {
        Vector3DI result;
        result.x = vector->x - op;
        result.y = vector->y - op;
        result.z = vector->z - op;
        return result;
    }

    Vector3DI Vector3DI_subtractVector(Vector3DI *vector, Vector3DI *op) {
        Vector3DI result;
        result.x = vector->x - op->x;
        result.y = vector->y - op->y;
        result.z = vector->z - op->z;
        return result;
    }

    //conflicting type with line 849

//     Vector3DI Vector3DI_multiplyInt(Vector3DI *vector, int op) {

//         Vector3DI result;

//         result.x = vector->x * (float)op;

//         result.y = vector->y * (float)op;

//         result.z = vector->z * (float)op;

//         return result;

//     }


    Vector3DI Vector3DI_multiplyFloat(Vector3DI *vector, float op) {
        Vector3DI result;
        result.x = vector->x * op;
        result.y = vector->y * op;
        result.z = vector->z * op;
        return result;

    }

    Vector3DI Vector3DI_multiplyVector(Vector3DI *vector, Vector3DI *op) {
        Vector3DI result;
        result.x = vector->x * op->x;
        result.y = vector->y * op->y;
        result.z = vector->z * op->z;
        return result;
    }

    Vector3DI Vector3DI_cross(Vector3DI *v1, Vector3DI *v2) {
        Vector3DI result;
        result.x = v1->y * v2->z - v1->z * v2->y;
        result.y = v1->z * v2->x - v1->x * v2->z;
        result.z = v1->x * v2->y - v1->y * v2->x;
        return result;
    }

    double Vector3DI_dot(Vector3DI *v1, Vector3DI *v2) {
        return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
    }
/*

    double Vector3DI_dist(Vector3DI *v1, Vector3DI *v2) {

        double dx = v1->x - v2->x;

        double dy = v1->y - v2->y;

        double dz = v1->z - v2->z;

        return sqrt(dx * dx + dy * dy + dz * dz);

    }*/


    double Vector3DI_distSq(Vector3DI *v1, Vector3DI *v2) {
        double dx = v1->x - v2->x;
        double dy = v1->y - v2->y;
        double dz = v1->z - v2->z;
        return dx * dx + dy * dy + dz * dz;
    }

/*

    Vector3DI* Vector3DI_Random(Vector3DI *v) {

    v->x = (float)rand() / RAND_MAX;

    v->y = (float)rand() / RAND_MAX;

    v->z = (float)rand() / RAND_MAX;

    return v;

}


    Vector3DI* Vector3DI_RandomWithRange(Vector3DI *v, Vector3DI a, Vector3DI b) {

        v->x = a.x + ((float)rand() / RAND_MAX) * (b.x - a.x);

        v->y = a.y + ((float)rand() / RAND_MAX) * (b.y - a.y);

        v->z = a.z + ((float)rand() / RAND_MAX) * (b.z - a.z);

        return v;

    }


    Vector3DI* Vector3DI_RandomWithBounds(Vector3DI *v, float x1, float x2, float y1, float y2, float z1, float z2) {

        v->x = x1 + ((float)rand() / RAND_MAX) * (x2 - x1);

        v->y = y1 + ((float)rand() / RAND_MAX) * (y2 - y1);

        v->z = z1 + ((float)rand() / RAND_MAX) * (z2 - z1);

        return v;

    }*/


    Vector3DI Vector3DI_RGBtoHSV(Vector3DI v) {
        // implementation not provided
        return v;
    }

    Vector3DI Vector3DI_HSVtoRGB(Vector3DI v) {
        // implementation not provided
        return v;
    }

    Vector3DI* Vector3DI_Normalize(Vector3DI *v) {
        // implementation not provided
        return v;
    }

    Vector3DI* Vector3DI_Clamp(Vector3DI *v, float a, float b) {
        // implementation not provided
        return v;
    }

    double Vector3DI_Length(Vector3DI v) {
        // implementation not provided
        return 0;
    }

	// Vector3DF Declarations

	#define VNAME		3DF
	#define VTYPE		float

	///////////Vector 3DF functions

	Vector3DF Vector3DF_init() {
		Vector3DF v;
		v.x = 0;
		v.y = 0;
		v.z = 0;
		return v;
	}

	void Vector3DF_destroy(Vector3DF *v) {
		// do nothing
	}

	void printVector3DF(const Vector3DF* vector) {
		printf("(%lf, %lf, %lf)", vector->x, vector->y, vector->z);
	}

	Vector3DF Vector3DF_init_with_values(float xa, float ya, float za) {
		Vector3DF v;
		v.x = xa;
		v.y = ya;
		v.z = za;
		return v;
	}

	void Set(Vector3DF* vec, float x, float y, float z) {
		vec->x = x;
		vec->y = y;
		vec->z = z;
	}

	Vector3DF Vector3DF_init_with_Vector2DC(Vector2DI *op) {
		Vector3DF v;
		v.x = op->x;
		v.y = op->y;
		v.z = 0;
		return v;
	}

	Vector3DF Vector3DF_init_with_Vector2DI(Vector2DI *op) {
		Vector3DF v;
		v.x = op->x;
		v.y = op->y;
		v.z = 0;
		return v;
	}

	Vector3DF Vector3DF_init_with_Vector2DF(Vector2DF *op) {
		Vector3DF v;
		v.x = op->x;
		v.y = op->y;
		v.z = 0;
		return v;
	}

	Vector3DF Vector3DF_init_with_Vector3DC(Vector3DC *op) {
		Vector3DF v;
		v.x = op->x;
		v.y = op->y;
		v.z = op->z;
		return v;
	}

	Vector3DF Vector3DF_init_with_Vector3DI(Vector3DI *op) {
		Vector3DF v;
		v.x = op->x;
		v.y = op->y;
		v.z = op->z;
		return v;
	}

	Vector3DF Vector3DF_init_with_Vector3DF(Vector3DF *op) {
		Vector3DF v;
		v.x = op->x;
		v.y = op->y;
		v.z = op->z;
		return v;
	}

	Vector3DF Vector3DF_init_with_Vector4DF(Vector4DF *op) {
		Vector3DF v;
		v.x = op->x;
		v.y = op->y;
		v.z = op->z;
		return v;
	}

	Vector3DF *Vector3DF_Set(Vector3DF *v, double xa, double ya, double za) {
		v->x = xa;
		v->y = ya;
		v->z = za;
		return v;
	}

	Vector3DF *Vector3DF_operator_equal_int(Vector3DF *v, int op) {
		v->x = op;
		v->y = op;
		v->z = op;
		return v;
	}

	Vector3DF *Vector3DF_operator_equal_double(Vector3DF *v, double op) {
		v->x = op;
		v->y = op;
		v->z = op;
		return v;
	}

	Vector3DF *Vector3DF_operator_equal_float(Vector3DF *v, float op) {
		v->x = op;
		v->y = op;
		v->z = op;
		return v;
	}

	Vector3DF *Vector3DF_operator_equal_Vector2DC(Vector3DF *v, Vector2DI *op) {
		v->x = op->x;
		v->y = op->y;
		v->z = 0;
		return v;
	}

	Vector3DF *Vector3DF_operator_equal_Vector2DI(Vector3DF *v, Vector2DI *op) {
		v->x = op->x;
		v->y = op->y;
		v->z = 0;
		return v;
	}

	Vector3DF *Vector3DF_operator_equal_Vector2DF(Vector3DF *v, Vector2DF *op) {
		v->x = op->x;
		v->y = op->y;
		v->z = 0;
		return v;
	}

	Vector3DF *Vector3DF_operator_equal_Vector3DC(Vector3DF *v, Vector3DC *op) {
		v->x = op->x;
		v->y = op->y;
		v->z = op->z;
		return v;
	}

	Vector3DF *Vector3DF_operator_equal_Vector3DI(Vector3DF *v, Vector3DI *op) {
		v->x = (float)op->x;
		v->y = (float)op->y;
		v->z = (float)op->z;
		return v;
	}

	Vector3DF *Vector3DF_operator_equal_Vector3DF(Vector3DF *v, Vector3DF *op) {
		v->x = op->x;
		v->y = op->y;
		v->z = op->z;
		return v;
	}

	Vector3DF *Vector3DF_operator_equal_Vector4DF(Vector3DF *v, Vector4DF *op) {
		v->x = op->x;
		v->y = op->y;
		v->z = op->z;
		return v;
	}

	Vector3DF Clamp(Vector3DF* vec, float a, float b) {
		Vector3DF result;
		result.x = (vec->x < a) ? a : ((vec->x > b) ? b : vec->x);
		result.y = (vec->y < a) ? a : ((vec->y > b) ? b : vec->y);
		result.z = (vec->z < a) ? a : ((vec->z > b) ? b : vec->z);
		return result;
	}
	Vector3DF *Vector3DF_addInt(Vector3DF *vector, int op) {
    vector->x += op;
    vector->y += op;
    vector->z += op;
    return vector;
	}

	Vector3DF *Vector3DF_addDouble(Vector3DF *vector, double op) {
		vector->x += op;
		vector->y += op;
		vector->z += op;
		return vector;
	}

	Vector3DF *Vector3DF_addVector2DC(Vector3DF *vector, Vector2DI *op) {
		vector->x += op->x;
		vector->y += op->y;
		return vector;
	}

	Vector3DF *Vector3DF_addVector2DI(Vector3DF *vector, Vector2DI *op) {
		vector->x += op->x;
		vector->y += op->y;
		return vector;
	}

	Vector3DF *Vector3DF_addVector2DF(Vector3DF *vector, Vector2DF *op) {
		vector->x += op->x;
		vector->y += op->y;
		return vector;
	}

	Vector3DF *Vector3DF_addVector3DC(Vector3DF *vector, Vector3DC *op) {
		vector->x += op->x;
		vector->y += op->y;
		vector->z += op->z;
		return vector;
	}

	Vector3DF *Vector3DF_addVector3DI(Vector3DF *vector, Vector3DI *op) {
		vector->x += op->x;
		vector->y += op->y;
		vector->z += op->z;
		return vector;
	}

	Vector3DF *Vector3DF_addVector3DF(Vector3DF *vector, Vector3DF *op) {
		vector->x += op->x;
		vector->y += op->y;
		vector->z += op->z;
		return vector;
	}

	Vector3DF *Vector3DF_addVector4DF(Vector3DF *vector, Vector4DF *op) {
		vector->x += op->x;
		vector->y += op->y;
		vector->z += op->z;
		return vector;
	}

	Vector3DF Vector3DF_subtractInt(Vector3DF *v, int op) {
		v->x -= op;
		v->y -= op;
		v->z -= op;
		return *v;
	}

	Vector3DF Vector3DF_subtractDouble(Vector3DF *v, double op) {
		v->x -= op;
		v->y -= op;
		v->z -= op;
		return *v;
	}

	Vector3DF Vector3DF_subtractVector2DC(Vector3DF *v, Vector2DI *op) {
		v->x -= op->x;
		v->y -= op->y;
		return *v;
	}

	Vector3DF Vector3DF_subtractVector2DI(Vector3DF *v, Vector2DI *op) {
		v->x -= op->x;
		v->y -= op->y;
		return *v;
	}

	Vector3DF Vector3DF_subtractVector2DF(Vector3DF *v, Vector2DF *op) {
		v->x -= op->x;
		v->y -= op->y;
		return *v;
	}

	Vector3DF Vector3DF_subtractVector3DC(Vector3DF *v, Vector3DC *op) {
		v->x -= op->x;
		v->y -= op->y;
		v->z -= op->z;
		return *v;
	}

	Vector3DF Vector3DF_subtractVector3DI(Vector3DF *v, Vector3DI *op) {
		v->x -= op->x;
		v->y -= op->y;
		v->z -= op->z;
		return *v;
	}

	Vector3DF Vector3DF_subtractVector3DF(Vector3DF *v, Vector3DF *op) {
		v->x -= op->x;
		v->y -= op->y;
		v->z -= op->z;
		return *v;
	}

	Vector3DF Vector3DF_subtractVector4DF(Vector3DF *v, Vector4DF *op) {
		v->x -= op->x;
		v->y -= op->y;
		v->z -= op->z;
		return *v;
	}

	Vector3DF *Vector3DF_multiplyInt(Vector3DF *v, int *op) {
		v->x *= *op;
		v->y *= *op;
		v->z *= *op;
		return v;
	}

	Vector3DF *Vector3DF_multiplyDouble(Vector3DF *v, double op) {
		v->x *= op;
		v->y *= op;
		v->z *= op;
		return v;
	}

	Vector3DF *Vector3DF_multiplyVector2DC(Vector3DF *v, Vector2DI *op) {
		v->x *= op->x;
		v->y *= op->y;
    return v;
	}

	Vector3DF *Vector3DF_multiplyVector2DI(Vector3DF *v, Vector2DI *op) {
    v->x *= op->x;
    v->y *= op->y;
    return v;
}

	Vector3DF *Vector3DF_multiplyVector2DF(Vector3DF *v, Vector2DF *op) {
		v->x *= op->x;
		v->y *= op->y;
		return v;
	}

	Vector3DF *Vector3DF_multiplyVector3DC(Vector3DF *v, Vector3DC *op) {
		v->x *= op->x;
		v->y *= op->y;
		v->z *= op->z;
		return v;
	}

	Vector3DF *Vector3DF_multiplyVector3DI(Vector3DF *v, Vector3DI *op) {
		v->x *= op->x;
		v->y *= op->y;
		v->z *= op->z;
		return v;
	}

	Vector3DF *Vector3DF_multiplyVector3DF(Vector3DF *v, Vector3DF *op) {
		v->x *= op->x;
		v->y *= op->y;
		v->z *= op->z;
		return v;
	}

	Vector3DF *Vector3DF_multiplyVector4DF(Vector3DF *v, Vector4DF *op) {
		v->x *= op->x;
		v->y *= op->y;
		v->z *= op->z;
		return v;
	}
/*
	Vector3DF *Vector3DF_multiplyMatrix4F(Vector3DF *v, Matrix4F *op) {
		// implementation MISSING
		return v;
	}

	Vector3DF *Vector3DF_multiplyMatrixF(Vector3DF *v, MatrixF *op) {
		// implementation MISSING
		return v;
	}
*/
	Vector3DF Vector3DF_divideByInt(Vector3DF *v, int op) {
		v->x /= op;
		v->y /= op;
		v->z /= op;
		return *v;
	}

	Vector3DF Vector3DF_divideByDouble(Vector3DF *v, double op) {
		v->x /= op;
		v->y /= op;
		v->z /= op;
		return *v;
	}

	Vector3DF Vector3DF_divideByVector2DC(Vector3DF *v, Vector2DC op) {
		v->x /= op.x;
		v->y /= op.y;
		return *v;
	}

	Vector3DF Vector3DF_divideByVector2DI(Vector3DF *v, Vector2DI op) {
		v->x /= op.x;
		v->y /= op.y;
		return *v;
	}

	Vector3DF Vector3DF_divideByVector2DF(Vector3DF *v, Vector2DF op) {
		v->x /= op.x;
		v->y /= op.y;
		return *v;
	}

	Vector3DF Vector3DF_divideByVector3DC(Vector3DF *v, Vector3DC op) {
		v->x /= op.x;
		v->y /= op.y;
		v->z /= op.z;
		return *v;
	}

	Vector3DF Vector3DF_divideByVector3DI(Vector3DF *v, Vector3DI op) {
		v->x /= op.x;
		v->y /= op.y;
		v->z /= op.z;
		return *v;
	}

	Vector3DF Vector3DF_divideByVector3DF(Vector3DF *v, Vector3DF op) {
		v->x /= op.x;
		v->y /= op.y;
		v->z /= op.z;
		return *v;
	}

	Vector3DF Vector3DF_divideByVector4DF(Vector3DF *v, Vector4DF op) {
		v->x /= op.x;
		v->y /= op.y;
		v->z /= op.z;
		return *v;
	}
/*
	Vector3DF Vector3DF_addInt(Vector3DF *vector, int op) {
		Vector3DF result;
		result.x = vector->x + (float)op;
		result.y = vector->y + (float)op;
		result.z = vector->z + (float)op;
		return result;
	}*/

	Vector3DF Vector3DF_addFloat(Vector3DF *vector, float op) {
		Vector3DF result;
		result.x = vector->x + op;
		result.y = vector->y + op;
		result.z = vector->z + op;
		return result;
	}

	Vector3DF Vector3DF_addVector(Vector3DF *vector, Vector3DF *op) {
		Vector3DF result;
		result.x = vector->x + op->x;
		result.y = vector->y + op->y;
		result.z = vector->z + op->z;
		return result;
	}

	//Redefinition of Vector_3DF_subtractInt???
	/*Vector3DF Vector3DF_subtractInt(Vector3DF *vector, int op) {
		Vector3DF result;
		result.x = vector->x - (float)op;
		result.y = vector->y - (float)op;
		result.z = vector->z - (float)op;
		return result;
	}*/

	Vector3DF Vector3DF_subtractFloat(Vector3DF *vector, float op) {
		Vector3DF result;
		result.x = vector->x - op;
		result.y = vector->y - op;
		result.z = vector->z - op;
		return result;
	}

	Vector3DF Vector3DF_subtractVector(Vector3DF *vector, Vector3DF *op) {
		Vector3DF result;
		result.x = vector->x - op->x;
		result.y = vector->y - op->y;
		result.z = vector->z - op->z;
		return result;
	}

	//conflicting type with line 849
// 	Vector3DF Vector3DF_multiplyInt(Vector3DF *vector, int op) {
// 		Vector3DF result;
// 		result.x = vector->x * (float)op;
// 		result.y = vector->y * (float)op;
// 		result.z = vector->z * (float)op;
// 		return result;
// 	}

	Vector3DF Vector3DF_multiplyFloat(Vector3DF *vector, float op) {
		Vector3DF result;
		result.x = vector->x * op;
		result.y = vector->y * op;
		result.z = vector->z * op;
		return result;
	}

	Vector3DF Vector3DF_multiplyVector(Vector3DF *vector, Vector3DF *op) {
		Vector3DF result;
		result.x = vector->x * op->x;
		result.y = vector->y * op->y;
		result.z = vector->z * op->z;
		return result;
	}

	Vector3DF Vector3DF_cross(Vector3DF *v1, Vector3DF *v2) {
		Vector3DF result;
		result.x = v1->y * v2->z - v1->z * v2->y;
		result.y = v1->z * v2->x - v1->x * v2->z;
		result.z = v1->x * v2->y - v1->y * v2->x;
		return result;
	}

	double Vector3DF_dot(Vector3DF *v1, Vector3DF *v2) {
		return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
	}
/*
	double Vector3DF_dist(Vector3DF *v1, Vector3DF *v2) {
		double dx = v1->x - v2->x;
		double dy = v1->y - v2->y;
		double dz = v1->z - v2->z;
		return sqrt(dx * dx + dy * dy + dz * dz);
	}*/

	double Vector3DF_distSq(Vector3DF *v1, Vector3DF *v2) {
		double dx = v1->x - v2->x;
		double dy = v1->y - v2->y;
		double dz = v1->z - v2->z;
		return dx * dx + dy * dy + dz * dz;
	}
/*
	Vector3DF* Vector3DF_Random(Vector3DF *v) {
    v->x = (float)rand() / RAND_MAX;
    v->y = (float)rand() / RAND_MAX;
    v->z = (float)rand() / RAND_MAX;
    return v;
}

	Vector3DF* Vector3DF_RandomWithRange(Vector3DF *v, Vector3DF a, Vector3DF b) {
		v->x = a.x + ((float)rand() / RAND_MAX) * (b.x - a.x);
		v->y = a.y + ((float)rand() / RAND_MAX) * (b.y - a.y);
		v->z = a.z + ((float)rand() / RAND_MAX) * (b.z - a.z);
		return v;
	}

	Vector3DF* Vector3DF_RandomWithBounds(Vector3DF *v, float x1, float x2, float y1, float y2, float z1, float z2) {
		v->x = x1 + ((float)rand() / RAND_MAX) * (x2 - x1);
		v->y = y1 + ((float)rand() / RAND_MAX) * (y2 - y1);
		v->z = z1 + ((float)rand() / RAND_MAX) * (z2 - z1);
		return v;
	}*/

	Vector3DF Vector3DF_RGBtoHSV(Vector3DF v) {
		// implementation not provided
		return v;
	}

	Vector3DF Vector3DF_HSVtoRGB(Vector3DF v) {
		// implementation not provided
		return v;
	}

	Vector3DF* Vector3DF_Normalize(Vector3DF *v) {
		// implementation not provided
		return v;
	}

	Vector3DF* Vector3DF_Clamp(Vector3DF *v, float a, float b) {
		// implementation not provided
		return v;
	}

	double Vector3DF_Length(Vector3DF v) {
		// implementation not provided
		return 0;
	}

/*
	VTYPE* Vector3DF_Data(Vector3DF *v) {
		float re = v->x;
		return *re;
	}*/
// 	/*/*class LUNA_CORE Vector3DF {
// 	public:
// 		VTYPE x, y, z;
//
//
// 		// Constructors/Destructors
// 		inline Vector3DF();
// 		inline ~Vector3DF();
// 		inline Vector3DF (const VTYPE xa, const VTYPE ya, const VTYPE za);
// 		inline Vector3DF (const Vector2DC &op);
// 		inline Vector3DF (const Vector2DI &op);
// 		inline Vector3DF (const Vector2DF &op);
// 		inline Vector3DF (const Vector3DC &op);
// 		inline Vector3DF (const Vector3DI &op);
// 		inline Vector3DF (const Vector3DF &op);
// 		inline Vector3DF (const Vector4DF &op);
//
// 		// Set Functions
// 		inline Vector3DF &Set (const double xa, const double ya, const double za);
//
// 		// Member Functions
// 		inline Vector3DF &operator= (const int op);
// 		inline Vector3DF &operator= (const double op);
// 		inline Vector3DF &operator= (const Vector2DC &op);
// 		inline Vector3DF &operator= (const Vector2DI &op);
// 		inline Vector3DF &operator= (const Vector2DF &op);
// 		inline Vector3DF &operator= (const Vector3DC &op);
// 		inline Vector3DF &operator= (const Vector3DI &op);
// 		inline Vector3DF &operator= (const Vector3DF &op);
// 		inline Vector3DF &operator= (const Vector4DF &op);*/*/
/*
		inline Vector3DF &operator+= (const int op);
		inline Vector3DF &operator+= (const double op);
		inline Vector3DF &operator+= (const Vector2DC &op);
		inline Vector3DF &operator+= (const Vector2DI &op);
		inline Vector3DF &operator+= (const Vector2DF &op);
		inline Vector3DF &operator+= (const Vector3DC &op);
		inline Vector3DF &operator+= (const Vector3DI &op);
		inline Vector3DF &operator+= (const Vector3DF &op);
		inline Vector3DF &operator+= (const Vector4DF &op);

		inline Vector3DF &operator-= (const int op);
		inline Vector3DF &operator-= (const double op);
		inline Vector3DF &operator-= (const Vector2DC &op);
		inline Vector3DF &operator-= (const Vector2DI &op);
		inline Vector3DF &operator-= (const Vector2DF &op);
		inline Vector3DF &operator-= (const Vector3DC &op);
		inline Vector3DF &operator-= (const Vector3DI &op);
		inline Vector3DF &operator-= (const Vector3DF &op);
		inline Vector3DF &operator-= (const Vector4DF &op);

		inline Vector3DF &operator*= (const int op);
		inline Vector3DF &operator*= (const double op);
		inline Vector3DF &operator*= (const Vector2DC &op);
		inline Vector3DF &operator*= (const Vector2DI &op);
		inline Vector3DF &operator*= (const Vector2DF &op);
		inline Vector3DF &operator*= (const Vector3DC &op);
		inline Vector3DF &operator*= (const Vector3DI &op);
		inline Vector3DF &operator*= (const Vector3DF &op);
		inline Vector3DF &operator*= (const Vector4DF &op);
		Vector3DF &operator*= (const Matrix4F &op);
		Vector3DF &operator*= (const MatrixF &op);				// see vector.cpp

		inline Vector3DF &operator/= (const int op);
		inline Vector3DF &operator/= (const double op);
		inline Vector3DF &operator/= (const Vector2DC &op);
		inline Vector3DF &operator/= (const Vector2DI &op);
		inline Vector3DF &operator/= (const Vector2DF &op);
		inline Vector3DF &operator/= (const Vector3DC &op);
		inline Vector3DF &operator/= (const Vector3DI &op);
		inline Vector3DF &operator/= (const Vector3DF &op);
		inline Vector3DF &operator/= (const Vector4DF &op);

		// Slow operations - require temporary variables
		inline Vector3DF operator+ (int op)			{ return Vector3DF(x+float(op), y+float(op), z+float(op)); }
		inline Vector3DF operator+ (float op)		{ return Vector3DF(x+op, y+op, z+op); }
		inline Vector3DF operator+ (Vector3DF &op)	{ return Vector3DF(x+op.x, y+op.y, z+op.z); }
		inline Vector3DF operator- (int op)			{ return Vector3DF(x-float(op), y-float(op), z-float(op)); }
		inline Vector3DF operator- (float op)		{ return Vector3DF(x-op, y-op, z-op); }
		inline Vector3DF operator- (Vector3DF &op)	{ return Vector3DF(x-op.x, y-op.y, z-op.z); }
		inline Vector3DF operator* (int op)			{ return Vector3DF(x*float(op), y*float(op), z*float(op)); }
		inline Vector3DF operator* (float op)		{ return Vector3DF(x*op, y*op, z*op); }
		inline Vector3DF operator* (Vector3DF &op)	{ return Vector3DF(x*op.x, y*op.y, z*op.z); }
		// --


		inline Vector3DF &Cross (const Vector3DC &v);
		inline Vector3DF &Cross (const Vector3DI &v);
		inline Vector3DF &Cross (const Vector3DF &v);

		inline double Dot(const Vector3DC &v);
		inline double Dot(const Vector3DI &v);
		inline double Dot(const Vector3DF &v);

		inline double Dist (const Vector2DC &v);
		inline double Dist (const Vector2DI &v);
		inline double Dist (const Vector2DF &v);
		inline double Dist (const Vector3DC &v);
		inline double Dist (const Vector3DI &v);
		inline double Dist (const Vector3DF &v);
		inline double Dist (const Vector4DF &v);

		inline double DistSq (const Vector2DC &v);
		inline double DistSq (const Vector2DI &v);
		inline double DistSq (const Vector2DF &v);
		inline double DistSq (const Vector3DC &v);
		inline double DistSq (const Vector3DI &v);
		inline double DistSq (const Vector3DF &v);
		inline double DistSq (const Vector4DF &v);

		inline Vector3DF &Random ()		{ x=float(rand())/RAND_MAX; y=float(rand())/RAND_MAX; z=float(rand())/RAND_MAX;  return *this;}
		inline Vector3DF &Random (Vector3DF a, Vector3DF b)		{ x=a.x+float(rand()*(b.x-a.x))/RAND_MAX; y=a.y+float(rand()*(b.y-a.y))/RAND_MAX; z=a.z+float(rand()*(b.z-a.z))/RAND_MAX;  return *this;}
		inline Vector3DF &Random (float x1,float x2, float y1, float y2, float z1, float z2)	{ x=x1+float(rand()*(x2-x1))/RAND_MAX; y=y1+float(rand()*(y2-y1))/RAND_MAX; z=z1+float(rand()*(z2-z1))/RAND_MAX;  return *this;}

		Vector3DF RGBtoHSV ();
		Vector3DF HSVtoRGB ();

		inline Vector3DF &Normalize (void);
        Vector3DF& Clamp (float a, float b);
		inline double Length (void);
		inline VTYPE *Data ();
	};*/

	#undef VNAME
	#undef VTYPE

// 	// Vector4DC Declarations
//
// 	#define VNAME		4DC
// 	#define VTYPE		unsigned char
//
// 	class LUNA_CORE Vector4DC {
// 	public:
// 		VTYPE x, y, z, w;
//
// 		inline Vector4DC &Set (const float xa, const float ya, const float za)	{ x = (VTYPE) xa; y= (VTYPE) ya; z=(VTYPE) za; w=1; return *this;}
// 		inline Vector4DC &Set (const float xa, const float ya, const float za, const float wa )	{ x =(VTYPE) xa; y= (VTYPE) ya; z=(VTYPE) za; w=(VTYPE) wa; return *this;}
// 		inline Vector4DC &Set (const VTYPE xa, const VTYPE ya, const VTYPE za)	{ x = (VTYPE) xa; y= (VTYPE) ya; z=(VTYPE) za; w=1; return *this;}
// 		inline Vector4DC &Set (const VTYPE xa, const VTYPE ya, const VTYPE za, const VTYPE wa )	{ x =(VTYPE) xa; y= (VTYPE) ya; z=(VTYPE) za; w=(VTYPE) wa; return *this;}
//
// 		// Constructors/Destructors
// 		inline Vector4DC();
// 		inline ~Vector4DC();
// 		inline Vector4DC (const VTYPE xa, const VTYPE ya, const VTYPE za, const VTYPE wa);
// 		inline Vector4DC (const Vector2DC &op);
// 		inline Vector4DC (const Vector2DI &op);
// 		inline Vector4DC (const Vector2DF &op);
// 		inline Vector4DC (const Vector3DC &op);
// 		inline Vector4DC (const Vector3DI &op);
// 		inline Vector4DC (const Vector3DF &op);
// 		inline Vector4DC (const Vector4DC &op);
// 		inline Vector4DC (const Vector4DF &op);
//
// 		// Member Functions
// 		inline Vector4DC &operator= ( const int op);
// 		inline Vector4DC &operator= ( const double op);
// 		inline Vector4DC &operator= ( const Vector2DC &op);
// 		inline Vector4DC &operator= ( const Vector2DI &op);
// 		inline Vector4DC &operator= ( const Vector2DF &op);
// 		inline Vector4DC &operator= ( const Vector3DC &op);
// 		inline Vector4DC &operator= ( const Vector3DI &op);
// 		inline Vector4DC &operator= ( const Vector3DF &op);
// 		inline Vector4DC &operator= ( const Vector4DC &op);
// 		inline Vector4DC &operator= ( const Vector4DF &op);
//
// 		inline Vector4DC &operator+= ( const int op);
// 		inline Vector4DC &operator+= ( const double op);
// 		inline Vector4DC &operator+= ( const Vector2DC &op);
// 		inline Vector4DC &operator+= ( const Vector2DI &op);
// 		inline Vector4DC &operator+= ( const Vector2DF &op);
// 		inline Vector4DC &operator+= ( const Vector3DC &op);
// 		inline Vector4DC &operator+= ( const Vector3DI &op);
// 		inline Vector4DC &operator+= ( const Vector3DF &op);
// 		inline Vector4DC &operator+= ( const Vector4DC &op);
// 		inline Vector4DC &operator+= ( const Vector4DF &op);
//
// 		inline Vector4DC &operator-= ( const int op);
// 		inline Vector4DC &operator-= ( const double op);
// 		inline Vector4DC &operator-= ( const Vector2DC &op);
// 		inline Vector4DC &operator-= ( const Vector2DI &op);
// 		inline Vector4DC &operator-= ( const Vector2DF &op);
// 		inline Vector4DC &operator-= ( const Vector3DC &op);
// 		inline Vector4DC &operator-= ( const Vector3DI &op);
// 		inline Vector4DC &operator-= ( const Vector3DF &op);
// 		inline Vector4DC &operator-= ( const Vector4DC &op);
// 		inline Vector4DC &operator-= ( const Vector4DF &op);
//
// 		inline Vector4DC &operator*= ( const int op);
// 		inline Vector4DC &operator*= ( const double op);
// 		inline Vector4DC &operator*= ( const Vector2DC &op);
// 		inline Vector4DC &operator*= ( const Vector2DI &op);
// 		inline Vector4DC &operator*= ( const Vector2DF &op);
// 		inline Vector4DC &operator*= ( const Vector3DC &op);
// 		inline Vector4DC &operator*= ( const Vector3DI &op);
// 		inline Vector4DC &operator*= ( const Vector3DF &op);
// 		inline Vector4DC &operator*= ( const Vector4DC &op);
// 		inline Vector4DC &operator*= ( const Vector4DF &op);
//
// 		inline Vector4DC &operator/= ( const int op);
// 		inline Vector4DC &operator/= ( const double op);
// 		inline Vector4DC &operator/= ( const Vector2DC &op);
// 		inline Vector4DC &operator/= ( const Vector2DI &op);
// 		inline Vector4DC &operator/= ( const Vector2DF &op);
// 		inline Vector4DC &operator/= ( const Vector3DC &op);
// 		inline Vector4DC &operator/= ( const Vector3DI &op);
// 		inline Vector4DC &operator/= ( const Vector3DF &op);
// 		inline Vector4DC &operator/= ( const Vector4DC &op);
// 		inline Vector4DC &operator/= ( const Vector4DF &op);
//
// 		// Slow operations - require temporary variables
// 		inline Vector4DC operator+ ( const int op);
// 		inline Vector4DC operator+ ( const float op);
// 		inline Vector4DC operator+ ( const Vector4DC &op);
// 		inline Vector4DC operator- ( const int op);
// 		inline Vector4DC operator- ( const float op);
// 		inline Vector4DC operator- ( const Vector4DC &op);
// 		inline Vector4DC operator* ( const int op);
// 		inline Vector4DC operator* ( const float op);
// 		inline Vector4DC operator* ( const Vector4DC &op);
// 		// --
//
// 		inline double Dot( const Vector4DF &v);
// 		inline double Dist ( const Vector4DF &v);
// 		inline double DistSq ( const Vector4DF &v);
// 		inline Vector4DC &Normalize (void);
// 		inline double Length (void);
//
// 		inline Vector4DC &Random ()		{ x=(VTYPE) float(rand()*255)/RAND_MAX; y=(VTYPE) float(rand()*255)/RAND_MAX; z=(VTYPE) float(rand()*255)/RAND_MAX; w = 1;  return *this;}
// 		inline VTYPE *Data (void);
// 	};
// 	#undef VNAME
// 	#undef VTYPE
//
//
// 	// Vector4DF Declarations
//
// 	#define VNAME		4DF
// 	#define VTYPE		float
//
// 	class LUNA_CORE Vector4DF {
// 	public:
// 		VTYPE x, y, z, w;
//
// 		inline Vector4DF &Set (const float xa, const float ya, const float za)	{ x =xa; y= ya; z=za; w=1; return *this;}
// 		inline Vector4DF &Set (const float xa, const float ya, const float za, const float wa )	{ x =xa; y= ya; z=za; w=wa; return *this;}
//
// 		// Constructors/Destructors
// 		inline Vector4DF();
// 		inline ~Vector4DF();
// 		inline Vector4DF (const VTYPE xa, const VTYPE ya, const VTYPE za, const VTYPE wa);
// 		inline Vector4DF (const Vector2DC &op);
// 		inline Vector4DF (const Vector2DI &op);
// 		inline Vector4DF (const Vector2DF &op);
// 		inline Vector4DF (const Vector3DC &op);
// 		inline Vector4DF (const Vector3DI &op);
// 		inline Vector4DF (const Vector3DF &op);
// 		inline Vector4DF (const Vector4DF &op);
//
// 		// Member Functions
// 		inline Vector4DF &operator= (const int op);
// 		inline Vector4DF &operator= (const double op);
// 		inline Vector4DF &operator= (const Vector2DC &op);
// 		inline Vector4DF &operator= (const Vector2DI &op);
// 		inline Vector4DF &operator= (const Vector2DF &op);
// 		inline Vector4DF &operator= (const Vector3DC &op);
// 		inline Vector4DF &operator= (const Vector3DI &op);
// 		inline Vector4DF &operator= (const Vector3DF &op);
// 		inline Vector4DF &operator= (const Vector4DF &op);
//
// 		inline Vector4DF &operator+= (const int op);
// 		inline Vector4DF &operator+= (const float op);
// 		inline Vector4DF &operator+= (const double op);
// 		inline Vector4DF &operator+= (const Vector2DC &op);
// 		inline Vector4DF &operator+= (const Vector2DI &op);
// 		inline Vector4DF &operator+= (const Vector2DF &op);
// 		inline Vector4DF &operator+= (const Vector3DC &op);
// 		inline Vector4DF &operator+= (const Vector3DI &op);
// 		inline Vector4DF &operator+= (const Vector3DF &op);
// 		inline Vector4DF &operator+= (const Vector4DF &op);
//
// 		inline Vector4DF &operator-= (const int op);
// 		inline Vector4DF &operator-= (const double op);
// 		inline Vector4DF &operator-= (const Vector2DC &op);
// 		inline Vector4DF &operator-= (const Vector2DI &op);
// 		inline Vector4DF &operator-= (const Vector2DF &op);
// 		inline Vector4DF &operator-= (const Vector3DC &op);
// 		inline Vector4DF &operator-= (const Vector3DI &op);
// 		inline Vector4DF &operator-= (const Vector3DF &op);
// 		inline Vector4DF &operator-= (const Vector4DF &op);
//
// 		inline Vector4DF &operator*= (const int op);
// 		inline Vector4DF &operator*= (const double op);
// 		inline Vector4DF &operator*= (const Vector2DC &op);
// 		inline Vector4DF &operator*= (const Vector2DI &op);
// 		inline Vector4DF &operator*= (const Vector2DF &op);
// 		inline Vector4DF &operator*= (const Vector3DC &op);
// 		inline Vector4DF &operator*= (const Vector3DI &op);
// 		inline Vector4DF &operator*= (const Vector3DF &op);
// 		inline Vector4DF &operator*= (const Vector4DF &op);
// 		Vector4DF &operator*= (const float* op );
// 		Vector4DF &operator*= (const Matrix4F &op);
// 		Vector4DF &operator*= (const MatrixF &op);				// see vector.cpp
//
// 		inline Vector4DF &operator/= (const int op);
// 		inline Vector4DF &operator/= (const double op);
// 		inline Vector4DF &operator/= (const Vector2DC &op);
// 		inline Vector4DF &operator/= (const Vector2DI &op);
// 		inline Vector4DF &operator/= (const Vector2DF &op);
// 		inline Vector4DF &operator/= (const Vector3DC &op);
// 		inline Vector4DF &operator/= (const Vector3DI &op);
// 		inline Vector4DF &operator/= (const Vector3DF &op);
// 		inline Vector4DF &operator/= (const Vector4DF &op);
//
// 		// Slow operations - require temporary variables
// 		inline Vector4DF operator+ (const int op)			{ return Vector4DF(x+float(op), y+float(op), z+float(op), w+float(op)); }
// 		inline Vector4DF operator+ (const float op)		{ return Vector4DF(x+op, y+op, z+op, w*op); }
// 		inline Vector4DF operator+ (const Vector4DF &op)	{ return Vector4DF(x+op.x, y+op.y, z+op.z, w+op.w); }
// 		inline Vector4DF operator- (const int op)			{ return Vector4DF(x-float(op), y-float(op), z-float(op), w-float(op)); }
// 		inline Vector4DF operator- (const float op)		{ return Vector4DF(x-op, y-op, z-op, w*op); }
// 		inline Vector4DF operator- (const Vector4DF &op)	{ return Vector4DF(x-op.x, y-op.y, z-op.z, w-op.w); }
// 		inline Vector4DF operator* (const int op)			{ return Vector4DF(x*float(op), y*float(op), z*float(op), w*float(op)); }
// 		inline Vector4DF operator* (const float op)		{ return Vector4DF(x*op, y*op, z*op, w*op); }
// 		inline Vector4DF operator* (const Vector4DF &op)	{ return Vector4DF(x*op.x, y*op.y, z*op.z, w*op.w); }
// 		// --
//
// 		inline Vector4DF &Set ( CLRVAL clr )	{
// 			x = RED(clr);		// (float( c      & 0xFF)/255.0)
// 			y = GRN(clr);		// (float((c>>8)  & 0xFF)/255.0)
// 			z = BLUE(clr);		// (float((c>>16) & 0xFF)/255.0)
// 			w = ALPH(clr);		// (float((c>>24) & 0xFF)/255.0)
// 			return *this;
// 		}
// 		inline Vector4DF& fromClr ( CLRVAL clr ) { return Set (clr); }
// 		inline CLRVAL toClr () { return (CLRVAL) COLORA( x, y, z, w ); }
//
// 		inline Vector4DF& Clamp ( float xc, float yc, float zc, float wc )
// 		{
// 			x = (x > xc) ? xc : x;
// 			y = (y > yc) ? yc : y;
// 			z = (z > zc) ? zc : z;
// 			w = (w > wc) ? wc : w;
// 			return *this;
// 		}
//
// 		inline Vector4DF &Cross (const Vector4DF &v);
//
// 		inline double Dot (const Vector4DF &v);
//
// 		inline double Dist (const Vector4DF &v);
//
// 		inline double DistSq (const Vector4DF &v);
//
// 		inline Vector4DF &Normalize (void);
// 		inline double Length (void);
//
// 		inline Vector4DF &Random ()		{ x=float(rand())/RAND_MAX; y=float(rand())/RAND_MAX; z=float(rand())/RAND_MAX; w = 1;  return *this;}
// 		inline VTYPE *Data (void);
// 	};

	#undef VNAME
	#undef VTYPE

    // Vector Code Definitions (Inlined)
	//#include "vector_inline.h"

#endif



