
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

        //#define		BUILD_CUDA			// CUDA - Visual Studio 2005 only (as of April 2009)

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
		double x, y;
	} Vector2DC;
	//class Vector2DI;
	typedef struct {
		double x, y;
	} Vector2DI;
	//class Vector2DF;
	typedef struct {
		double x, y;
	} Vector2DF;
	//class Vector3DC;
	typedef struct {
		double x, y, z;
	} Vector3DC;
	//class Vector3DI;
	typedef struct {
		double x, y, z;
	} Vector3DI;
	//class Vector3DF;
	typedef struct {
		double x, y, z;
	} Vector3DF;
	//class Vector4DC;
	typedef struct {
		double x, y, z;
	} Vector4DC;
	//class Vector4DF;
	typedef struct {
		double x, y, z, w;
	} Vector4DF;
	//class MatrixF;
	//class Matrix4F;


		// Constructors/Destructors
		inline Vector3DI();
		inline ~Vector3DI();
		inline Vector3DI (const VTYPE xa, const VTYPE ya, const VTYPE za);
		inline Vector3DI (const Vector2DC &op);
		inline Vector3DI (const Vector2DI &op);
		inline Vector3DI (const Vector2DF &op);
		inline Vector3DI (const Vector3DC &op);
		inline Vector3DI (const Vector3DI &op);
		inline Vector3DI (const Vector3DF &op);
		inline Vector3DI (const Vector4DF &op);

		// Set Functions
		inline Vector3DI &Set (const int xa, const int ya, const int za);

		// Member Functions
		inline Vector3DI &operator= (const Vector2DC &op);
		inline Vector3DI &operator= (const Vector2DI &op);
		inline Vector3DI &operator= (const Vector2DF &op);
		inline Vector3DI &operator= (const Vector3DC &op);
		inline Vector3DI &operator= (const Vector3DI &op);
		inline Vector3DI &operator= (const Vector3DF &op);
		inline Vector3DI &operator= (const Vector4DF &op);

		inline Vector3DI &operator+= (const Vector2DC &op);
		inline Vector3DI &operator+= (const Vector2DI &op);
		inline Vector3DI &operator+= (const Vector2DF &op);
		inline Vector3DI &operator+= (const Vector3DC &op);
		inline Vector3DI &operator+= (const Vector3DI &op);
		inline Vector3DI &operator+= (const Vector3DF &op);
		inline Vector3DI &operator+= (const Vector4DF &op);

		inline Vector3DI &operator-= (const Vector2DC &op);
		inline Vector3DI &operator-= (const Vector2DI &op);
		inline Vector3DI &operator-= (const Vector2DF &op);
		inline Vector3DI &operator-= (const Vector3DC &op);
		inline Vector3DI &operator-= (const Vector3DI &op);
		inline Vector3DI &operator-= (const Vector3DF &op);
		inline Vector3DI &operator-= (const Vector4DF &op);

		inline Vector3DI &operator*= (const Vector2DC &op);
		inline Vector3DI &operator*= (const Vector2DI &op);
		inline Vector3DI &operator*= (const Vector2DF &op);
		inline Vector3DI &operator*= (const Vector3DC &op);
		inline Vector3DI &operator*= (const Vector3DI &op);
		inline Vector3DI &operator*= (const Vector3DF &op);
		inline Vector3DI &operator*= (const Vector4DF &op);

		inline Vector3DI &operator/= (const Vector2DC &op);
		inline Vector3DI &operator/= (const Vector2DI &op);
		inline Vector3DI &operator/= (const Vector2DF &op);
		inline Vector3DI &operator/= (const Vector3DC &op);
		inline Vector3DI &operator/= (const Vector3DI &op);
		inline Vector3DI &operator/= (const Vector3DF &op);
		inline Vector3DI &operator/= (const Vector4DF &op);

		inline Vector3DI operator+ (const int op)			{ return Vector3DI(x+(VTYPE) op, y+(VTYPE) op, z+(VTYPE) op); }
		inline Vector3DI operator+ (const float op)		{ return Vector3DI(x+(VTYPE) op, y+(VTYPE) op, z+(VTYPE) op); }
		inline Vector3DI operator+ (const Vector3DI &op)	{ return Vector3DI(x+op.x, y+op.y, z+op.z); }
		inline Vector3DI operator- (const int op)			{ return Vector3DI(x-(VTYPE) op, y-(VTYPE) op, z-(VTYPE) op); }
		inline Vector3DI operator- (const float op)		{ return Vector3DI(x-(VTYPE) op, y-(VTYPE) op, z-(VTYPE) op); }
		inline Vector3DI operator- (const Vector3DI &op)	{ return Vector3DI(x-op.x, y-op.y, z-op.z); }
		inline Vector3DI operator* (const int op)			{ return Vector3DI(x*(VTYPE) op, y*(VTYPE) op, z*(VTYPE) op); }
		inline Vector3DI operator* (const float op)		{ return Vector3DI(x*(VTYPE) op, y*(VTYPE) op, z*(VTYPE) op); }
		inline Vector3DI operator* (const Vector3DI &op)	{ return Vector3DI(x*op.x, y*op.y, z*op.z); }

		inline Vector3DI &Cross (const Vector3DC &v);
		inline Vector3DI &Cross (const Vector3DI &v);
		inline Vector3DI &Cross (const Vector3DF &v);

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

		inline Vector3DI &Normalize (void);
		inline double Length (void);
		inline VTYPE *Data (void);
	};*/

	#undef VNAME
	#undef VTYPE

#endif



