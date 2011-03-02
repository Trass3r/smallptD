/**
 * smallptD
 * a tiny path tracer in D
 * (C) 2010-11 Trass3r
 *
 * originally written in C++ by Kevin Beason, 2008
 * Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
 * Remove "-fopenmp" for g++ version < 4.2
 * Usage: time ./smallpt 5000 && xv image.ppm
 */
module smallpt;

import std.conv;
import std.math;
import std.random;
import std.stdio;

struct Vec(T)
{
	T x=0, y=0, z=0;									// position, also color (r,g,b)
	Vec opAdd(const Vec b) const { return Vec(x+b.x,y+b.y,z+b.z); }
	Vec opSub(const Vec b) const { return Vec(x-b.x,y-b.y,z-b.z); }
	Vec opMul(T b) const { return Vec(x*b,y*b,z*b); }
	Vec opMul(const Vec b) const { return Vec(x*b.x,y*b.y,z*b.z); }
	ref Vec norm(){ this = this * (1/sqrt(x*x+y*y+z*z)); return this;}
	T dot(const Vec b) const { return x*b.x+y*b.y+z*b.z; } // cross:
	Vec opMod(const Vec b) const {return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
}
alias Vec!double Vec3;

unittest
{
	Vec3 zero, one = Vec3(1,1,1), v1=Vec3(-1,0,2), v2=Vec3(10,-20,5.5);
	assert(zero+one == one);
	assert(zero*one == zero);
	assert(v1+v2 == Vec3(9,-20,7.5));
	assert(v1-v2 == Vec3(-11,20,-3.5));
	assert(v1.dot(v2) == 1);
	assert(v1*v2 == Vec3(-10,-0.,11));
	v1.norm();
	auto t = v1 - Vec3(-0.44721359549995793928183473374626, 0, 0.89442719099991587856366946749251);
	assert(sqrt(t.x*t.x + t.y*t.y + t.z*t.z) <= 0.001);
}
struct Ray
{
	Vec3 o, d;
}
enum Refl { DIFF, SPEC, REFR };	// material types, used in radiance()
struct Sphere
{
	double rad=0;			 // radius
	Vec3 p, e, c;			// position, emission, color
	Refl refl;			// reflection type (DIFFuse, SPECular, REFRactive)
	double intersect(const ref Ray r) const
	{ // returns distance, 0 if nohit
		Vec3 op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		double t=0, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
		if (det<0)
			return 0;
		else
			det=sqrt(det);
		return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
	}
}
__gshared immutable Sphere[9] spheres = [ //Scene: radius, position, emission, color, material
	Sphere(1e5, Vec3( 1e5+1,40.8,81.6),	Vec3(),			Vec3(.75,.25,.25),	Refl.DIFF),//Left
	Sphere(1e5, Vec3(-1e5+99,40.8,81.6),	Vec3(),			Vec3(.25,.25,.75),	Refl.DIFF),//Rght
	Sphere(1e5, Vec3(50,40.8, 1e5),		Vec3(),			Vec3(.75,.75,.75),	Refl.DIFF),//Back
	Sphere(1e5, Vec3(50,40.8,-1e5+170),	Vec3(),			Vec3(),				Refl.DIFF),//Frnt
	Sphere(1e5, Vec3(50, 1e5, 81.6),		Vec3(),			Vec3(.75,.75,.75),	Refl.DIFF),//Botm
	Sphere(1e5, Vec3(50,-1e5+81.6,81.6),	Vec3(),			Vec3(.75,.75,.75),	Refl.DIFF),//Top
	Sphere(16.5,Vec3(27,16.5,47),		Vec3(),			Vec3(1,1,1)*.999,	Refl.SPEC),//Mirr
	Sphere(16.5,Vec3(73,16.5,78),		Vec3(),			Vec3(1,1,1)*.999,	Refl.REFR),//Glas
	Sphere(600, Vec3(50,681.6-.27,81.6),	Vec3(12,12,12),	Vec3(),				Refl.DIFF) //Lite
];
T clamp(T)(T x){ return x<0 ? 0 : x>1 ? 1 : x; }
int toInt(double x){ return cast(int)(pow(clamp(x),1/2.2)*255+.5); }
bool intersect(const ref Ray r, out double t, out size_t id)
{
	double d=0, inf=t=1e20;
	foreach_reverse(i, sphere; spheres)
	{
		d=sphere.intersect(r);
		if(d > 0 && d<t)
			{t=d; id=i;}
	}
	return t<inf;
}
Vec radiance(const ref Ray r, int depth, ushort *Xi = null)
{
	double t=0;	// distance to intersection
	size_t id=0;	// id of intersected object
	if (!intersect(r, t, id))
		return Vec3(); // if miss, return black
	const Sphere obj = spheres[id];				// the hit object
	Vec3 x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1, f=obj.c;
	double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
	if (++depth>5)
		if (uniform(0.0, 1.0)<p)
			f=f*(1/p);
		else
			return obj.e; //R.R.

	if (obj.refl == Refl.DIFF)
	{									// Ideal DIFFUSE reflection
		double r1=2*PI*uniform(0.0, 1.0), r2=uniform(0.0, 1.0), r2s=sqrt(r2);
		Vec3 w=nl, u=((fabs(w.x)>.1?Vec3(0,1):Vec3(1))%w).norm(), v=w%u;
		Vec3 d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
		return obj.e + f * radiance(Ray(x,d),depth,Xi);
	}
	else if (obj.refl == Refl.SPEC)						// Ideal SPECULAR reflection
		return obj.e + f * radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi);

	Ray reflRay = Ray(x, r.d-n*2*n.dot(r.d));		 // Ideal dielectric REFRACTION
	bool into = n.dot(nl)>0;								// Ray from outside going in?
	double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t=0;

	if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)		// Total internal reflection
		return obj.e + f * radiance(reflRay,depth,Xi);

	Vec3 tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
	double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
	double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
	return obj.e + f * (depth>2 ? (uniform(0.0, 1.0)<P ?	 // Russian roulette
		radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
		radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}
void main(string[] args)
{
	const int w=1024, h=768, samps = args.length==2 ? to!int(args[1])/4 : 1; // # samples
	Ray cam = Ray(Vec3(50,52,295.6), Vec3(0,-0.042612,-1).norm()); // cam pos, dir
	Vec3 cx = Vec3(w*.5135/h), cy=(cx%cam.d).norm()*.5135, r;
	Vec3[] c = new Vec3[w*h];
//pragma omp parallel for schedule(dynamic, 1) private(r)			 // OpenMP
	for (int y=0; y<h; y++)
	{											 // Loop over image rows
		stderr.writef("\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
//		ushort[3] Xi = [0,0,y*y*y];
		for (ushort x=0; x<w; x++)	 // Loop cols
			for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)		 // 2x2 subpixel rows
				for (int sx=0; sx<2; sx++, r=Vec3())
				{				// 2x2 subpixel cols
					for (int s=0; s<samps; s++)
					{
						double r1=2*uniform(0.0, 1.0), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
						double r2=2*uniform(0.0, 1.0), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
						Vec3 d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
										cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
						r = r + radiance(Ray(cam.o+d*140,d.norm()),0)*(1./samps);
					} // Camera rays are pushed ^^^^^ forward to start in interior
					c[i] = c[i] + Vec3(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
				}
	}
	File f = File("image.ppm", "w");				 // Write image to PPM file.
	f.writef("P3\n%d %d\n%d\n", w, h, 255);
	for (int i=0; i<w*h; i++)
		f.writef("%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	f.close();
}