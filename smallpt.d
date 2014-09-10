/**
 * smallptD
 * a tiny path tracer in D
 * (C) 2010-11 Trass3r
 *
 * originally written in C++ by Kevin Beason, 2008
 * Make : ldc2 -m64 -O3 -release -vectorize-slp-aggressive -g smallpt.d
 * Usage: time ./smallpt 5000 && xv image.ppm
 */
module smallpt;

import core.stdc.tgmath;
import std.random;

enum float PI = 0x1.921fb54442d18469898cc51701b84p+1L;

align(16) struct Vec(T)
{
	T x=0, y=0, z=0, w=0;			// position, also color (r,g,b)
	invariant() { assert(w == 0); }
	Vec opBinary(string op)(const Vec b) const { return mixin("Vec(x"~op~"b.x,y"~op~"b.y,z"~op~"b.z,w"~op~"b.w)"); }
	Vec opBinary(string op:"*")(T b) const { return Vec(x*b,y*b,z*b); }
	ref Vec norm(){ this = this * (1/sqrt(x*x+y*y+z*z+w*w)); return this;}
	T dot(const Vec b) const { return x*b.x+y*b.y+z*b.z+w*b.w; } // cross:
	Vec opMod(const Vec b) const {return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
}
alias Vec!float Vec3;

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
	assert(sqrt(t.x*t.x + t.y*t.y + t.z*t.z) <= 0.001f);
}
struct Ray
{
	Vec3 o, d;
}
enum Refl { DIFF, SPEC, REFR };	// material types, used in radiance()
struct Sphere
{
	Vec3 p, e, c;			// position, emission, color
	float rad=0;			 // radius
	Refl refl;			// reflection type (DIFFuse, SPECular, REFRactive)
	float intersect(const ref Ray r) const
	{ // returns distance, 0 if nohit
		Vec3 op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		float t=0, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
		if (det<0)
			return 0;
		else
			det=sqrt(det);
		return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
	}
}
__gshared immutable Sphere[9] spheres = [ //Scene: position, emission, color, radius, material
	Sphere(Vec3( 1e5+1,40.8,81.6), Vec3(),         Vec3(0.75,0.25,0.25), 1e5, Refl.DIFF),//Left
	Sphere(Vec3(-1e5+99,40.8,81.6),Vec3(),         Vec3(0.25,0.25,0.75), 1e5, Refl.DIFF),//Rght
	Sphere(Vec3(50,40.8, 1e5),     Vec3(),         Vec3(0.75,0.75,0.75), 1e5, Refl.DIFF),//Back
	Sphere(Vec3(50,40.8,-1e5+170), Vec3(),         Vec3(),               1e5, Refl.DIFF),//Frnt
	Sphere(Vec3(50, 1e5, 81.6),    Vec3(),         Vec3(0.75,0.75,0.75), 1e5, Refl.DIFF),//Botm
	Sphere(Vec3(50,-1e5+81.6,81.6),Vec3(),         Vec3(0.75,0.75,0.75), 1e5, Refl.DIFF),//Top
	Sphere(Vec3(27,16.5,47),       Vec3(),         Vec3(1,1,1)*0.999f,  16.5, Refl.SPEC),//Mirr
	Sphere(Vec3(73,16.5,78),       Vec3(),         Vec3(1,1,1)*0.999f,  16.5, Refl.REFR),//Glas
	Sphere(Vec3(50,81.6-16.5,81.6),Vec3(4,4,4)*100,Vec3(),               1.5, Refl.DIFF),//Lite
];
T clamp(T)(T x){ return x<0 ? 0 : (x>1 ? 1 : x); }
int toInt(float x){ return cast(int)(pow(clamp(x),1/2.2f)*255+0.5f); }
float fabs(float x) { return x<0 ? -x : x; }
bool intersect(const Ray r, out float t, out size_t id)
{
	float inf = t = 1e20;
	foreach(i, sphere; spheres)
	{
		float d = sphere.intersect(r);
		if(d > 0 && d<t)
			{t=d; id=i;}
	}
	return t<inf;
}
Vec3 radiance(const Ray r, int depth, ushort *Xi = null, int E=1)
{
	float t=0;	// distance to intersection
	size_t id=0;	// id of intersected object
	if (!intersect(r, t, id))
		return Vec3(); // if miss, return black
	const Sphere obj = spheres[id];				// the hit object
	Vec3 x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1, f=obj.c;
	float p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
	if (++depth>5 || !p)
		if (uniform(0.0f, 1.0f)<p)
			f=f*(1/p);
		else
			return obj.e * E;

	if (obj.refl == Refl.DIFF)
	{									// Ideal DIFFUSE reflection
		float r1=2*PI*uniform(0.0f, 1.0f), r2=uniform(0.0f, 1.0f), r2s=sqrt(r2);
		Vec3 w=nl, u=((fabs(w.x)>0.1f?Vec3(0,1):Vec3(1))%w).norm(), v=w%u;
		Vec3 d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();

		// Loop over any lights
		Vec3 e;
		foreach (i, s; spheres)
		{
			if (s.e.x<=0 && s.e.y<=0 && s.e.z<=0)
				continue; // skip non-lights

			Vec3 sw=s.p-x, su=((fabs(sw.x)>.1?Vec3(0,1):Vec3(1))%sw).norm(), sv=sw%su;
			float cos_a_max = sqrt(1-s.rad*s.rad/(x-s.p).dot(x-s.p));
			float eps1 = uniform(0.0f, 1.0f), eps2 = uniform(0.0f, 1.0f);
			float cos_a = 1-eps1+eps1*cos_a_max;
			float sin_a = sqrt(1-cos_a*cos_a);
			float phi = 2*PI*eps2;
			Vec3 l = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;
			l.norm();
			if (intersect(Ray(x,l), t, id) && id==i){	// shadow ray
				float omega = 2*PI*(1-cos_a_max);
				e = e + f * (s.e*l.dot(nl)*omega)*(1/PI);  // 1/pi for brdf
			}
		}

		return obj.e*E+e+f * (radiance(Ray(x,d),depth,Xi,0));
	}
	else if (obj.refl == Refl.SPEC)						// Ideal SPECULAR reflection
		return obj.e + f * radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi);

	Ray reflRay = Ray(x, r.d-n*2*n.dot(r.d));		 // Ideal dielectric REFRACTION
	bool into = n.dot(nl)>0;								// Ray from outside going in?
	float nc=1, nt=1.5f, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t=0;

	if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)		// Total internal reflection
		return obj.e + f * radiance(reflRay,depth,Xi);

	Vec3 tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
	float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
	float Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=0.25f+0.5f*Re,RP=Re/P,TP=Tr/(1-P);
	return obj.e + f * (depth>2 ? (uniform(0.0f, 1.0f)<P ?	 // Russian roulette
		radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
		radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}
void main(string[] args)
{
	import std.conv;
	import std.stdio;

	const int w=1024, h=768, samps = args.length==2 ? to!int(args[1])/4 : 1; // # samples
	Ray cam = Ray(Vec3(50,52,295.6f), Vec3(0,-0.042612f,-1).norm()); // cam pos, dir
	Vec3 cx = Vec3(w*0.5135f/h), cy=(cx%cam.d).norm()*0.5135f, r;
	Vec3[] c = new Vec3[w*h];
//pragma omp parallel for schedule(dynamic, 1) private(r)			 // OpenMP
	for (int y=0; y<h; y++)
	{											 // Loop over image rows
		stderr.writef("\rRendering (%d spp) %5.2f%%",samps*4,100.0f*y/(h-1));
//		ushort[3] Xi = [0,0,y*y*y];
		for (ushort x=0; x<w; ++x)	 // Loop cols
			for (int sy=0, i=(h-y-1)*w+x; sy<2; ++sy)		 // 2x2 subpixel rows
				for (int sx=0; sx<2; ++sx, r=Vec3())
				{				// 2x2 subpixel cols
					for (int s=0; s<samps; s++)
					{
						float r1=2*uniform(0.0f, 1.0f), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
						float r2=2*uniform(0.0f, 1.0f), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
						Vec3 d = cx*( ( (sx+0.5f + dx)/2 + x)/w - 0.5f) +
										cy*( ( (sy+0.5f + dy)/2 + y)/h - 0.5f) + cam.d;
						r = r + radiance(Ray(cam.o+d*140,d.norm()),0)*(1.0f/samps);
					} // Camera rays are pushed ^^^^^ forward to start in interior
					c[i] = c[i] + Vec3(clamp(r.x),clamp(r.y),clamp(r.z))*0.25f;
				}
	}
	File f = File("image.ppm", "w");				 // Write image to PPM file.
	f.writef("P3\n%d %d\n%d\n", w, h, 255);
	for (int i=0; i<w*h; i++)
		f.writef("%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	f.close();
}