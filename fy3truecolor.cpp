// fy3truecolor.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <string>
#include <iostream>
#include <vector>
#include "gdal_priv.h"
#include <cmath>
#include "../../wftools.h"


using namespace std ;
 
#define DEG2RAD	0.0174532925199		/* PI/180 */
#define UO3	0.319
#define UH2O	2.93
#define REFLMIN -0.01
#define REFLMAX  1.6
 
#define MAXSOLZ 86.5
#define MAXAIRMASS 18
#define	SCALEHEIGHT 8000
#define MAXNUMSPHALBVALUES	3000
#define	TAUMAX	0.3
#define	TAUSTEP4SPHALB	(TAUMAX / (float)MAXNUMSPHALBVALUES)


string gdaltranslate = "gdal_translate";


/// demarr is a 2darry of resolution 0.083 degree 4320x2160
short getDemByLonLat( short* demarr , int demxsize, int demysize , float lon , float lat)
{
	int ix = (int)(lon*12+2160+0.5) ;
	int iy = (int)(lat*(-12)+1080+0.5) ;
	return demarr[iy*demxsize+ix] ;
}




double fintexp1(float tau)
{
	double xx, xftau;
	int i;
	const float a[6] = {-.57721566, 0.99999193,-0.24991055,
					0.05519968,-0.00976004, 0.00107857};
	xx = a[0];
	xftau = 1.;
	for (i=1; i<6; i++) {
		xftau *= tau;
		xx += a[i] * xftau;
	}
	return xx - logf(tau);
}

//! to compute exponential integrads for argument tau. see wiki for detail.
double fintexp3(float tau)
{
	return (expf(-tau) * (1. - tau) + tau * tau * fintexp1(tau)) / 2.;
}

//! to compute the spherical albedo of molecular layer. see 6s manual part2 page96.
/// 计算大气球反照率
/// @param tau Rayleigh光学厚度
/// @return 返回大气球反照率
float csalbr(float tau)
{
  return (3 * tau - fintexp3(tau) * (4 + 2 * tau) + 2 * expf(-tau)) / (4 + 3 * tau);
}

/// 计算Rayleigh散射反射率
/// @param phi 太阳方位角与卫星方位角的差值，单位角度[-180~+180]
/// @param muv 卫星天顶角Sensor Zenith Angle的余弦值,cos(SensorZenithAngle)
/// @param mus 太阳天顶角Solar Zenith Angle 的余弦值,cos(SolarZenithAngle)
/// @param taur 输入Rayleigh光学厚度，taur[0]为band1（蓝），taur[1]为band2（绿），taur[2]为band3（红）
/// @param rhoray 返回三个波段Rayleigh散射反射率，对应波段同上
/// @param trup 返回三个波段Rayleigh上行透过率，对应波段同上
/// @param trdown 返回三个波段Rayleigh下行透过率，对应波段同上
void chand(
	float phi, 
	float muv, 
	float mus, 
	float *taur, 
	float *rhoray, 
	float *trup, 
	float *trdown
	)
{
	const double xfd=0.958725775;
	const float xbeta2=0.5;
	float pl[5];
	double fs01,fs02,fs0, fs1,fs2;
	const float as0[10] = {0.33243832, 0.16285370, -0.30924818, -0.10324388, 0.11493334,
						   -6.777104e-02, 1.577425e-03, -1.240906e-02, 3.241678e-02, -3.503695e-02};
	const float as1[2] = {.19666292, -5.439061e-02};
	const float as2[2] = {.14545937,-2.910845e-02};
	float phios,xcosf1,xcosf2,xcosf3;
	float xph1,xph2,xph3,xitm1,xitm2;
	float xlntaur,xitot1,xitot2,xitot3;
	int i,ib;

	phios = phi + 180;
	xcosf1 = 1.;
	xcosf2 = cosf(phios * DEG2RAD);
	xcosf3 = cosf(2 * phios * DEG2RAD);
	xph1 = 1 + (3 * mus * mus - 1) * (3 * muv * muv - 1) * xfd / 8.;
	xph2 = - xfd * xbeta2 * 1.5 * mus * muv * sqrtf(1 - mus * mus) * sqrtf(1 - muv * muv);
	xph3 =   xfd * xbeta2 * 0.375 * (1 - mus * mus) * (1 - muv * muv);
	pl[0] = 1.;
	pl[1] = mus + muv;
	pl[2] = mus * muv;
	pl[3] = mus * mus + muv * muv;
	pl[4] = mus * mus * muv * muv;
	fs01 = fs02 = 0;
	for (i=0; i<5; i++) fs01 += pl[i] * as0[i];
	for (i=0; i<5; i++) fs02 += pl[i] * as0[5 + i];
	for (ib=0; ib<3; ib++) {
		xlntaur = logf(taur[ib]);
		fs0 = fs01 + fs02 * xlntaur;
		fs1 = as1[0] + xlntaur * as1[1];
		fs2 = as2[0] + xlntaur * as2[1];
		trdown[ib] = expf(-taur[ib]/mus);
		trup[ib]   = expf(-taur[ib]/muv);
		xitm1 = (1 - trdown[ib] * trup[ib]) / 4. / (mus + muv);
		xitm2 = (1 - trdown[ib]) * (1 - trup[ib]);
		xitot1 = xph1 * (xitm1 + xitm2 * fs0);
		xitot2 = xph2 * (xitm1 + xitm2 * fs1);
		xitot3 = xph3 * (xitm1 + xitm2 * fs2);
		rhoray[ib] = xitot1 * xcosf1 + xitot2 * xcosf2 * 2 + xitot3 * xcosf3 * 2;
	}
}

/// 计算FY3D前三个可见光波段的rayleigh散射反射率和大气透过率
/// @param mus 太阳天顶角Solar Zenith Angle 的余弦值
/// @param muv 卫星天顶角Sensor Zenith Angle的余弦值
/// @param phi 太阳方位角与卫星方位角的差值
/// @param height 像元的海拔高度（米），不能小于0，低于海平面以0表示。
/// @param sphalb 返回三个波段大气球反照率 sphalb[0]为band1（蓝） sphalb[1]为band2（绿） sphalb[2]为band3（红）
/// @param rhoray 返回三个波段Rayleigh散射反射率，对应波段同上
/// @param TtotraytH2O 返回三个波段的Rayleigh和水汽综合透过率，对应波段同上
/// @param tOG 返回三个波段的O2和O3综合透过率，对应波段同上
/// @return 返回值为0，暂时没有用到
int getatmvariables(
	float mus, 
	float muv, 
	float phi, 
	short height, 
	float *sphalb, 
	float *rhoray, 
	float *TtotraytH2O, 
	float *tOG
	)
{
	double m, Ttotrayu, Ttotrayd, tO3, tO2, tH2O, psurfratio;
	int j, ib;
	///modis const float aH2O[Nbands]={ -5.60723, -5.25251, 0, 0, -6.29824, -7.70944, -3.91877 };
	///modis const float bH2O[Nbands]={ 0.820175, 0.725159, 0, 0, 0.865732, 0.966947, 0.745342 };
	///modis const float aO3[Nbands]={ 0.0715289, 0, 0.00743232, 0.089691, 0, 0, 0 };
	///modis const float taur0[Nbands] = { 0.05100, 0.01631, 0.19325, 0.09536, 0.00366, 0.00123, 0.00043 };

	/// FY3B band1(mod3) band2(mod4) band3(mod1)
	const float aH2O[3] = { 0 , 0 , -5.60723 } ;
	const float bH2O[3] = { 0 , 0 ,  0.820175 } ;
	const float aO3[3]={ 0.00743232, 0.089691 , 0.0715289  };
	const float taur0[3] = { 0.19325, 0.09536, 0.05100 };

	float taur[3], trup[3], trdown[3];
	static float sphalb0[MAXNUMSPHALBVALUES];
	static bool first_time=true;

	if (first_time) {
		sphalb0[0] = 0;
		for(j=1; j<MAXNUMSPHALBVALUES; j++)		/* taur <= 0.3 for bands 1 to 7 (including safety margin for height<~0) */
			sphalb0[j] = csalbr(j * TAUSTEP4SPHALB);//! compute molcular spherical albedo from tau 0.0 to 0.3.
		first_time = false;
	}

	m = 1 / mus + 1 / muv;
	if (m > MAXAIRMASS) return -1;
	psurfratio = expf(-height / (float)SCALEHEIGHT);
	for (ib=0; ib<3; ib++)
		taur[ib] = taur0[ib] * psurfratio;

	chand(phi, muv, mus, taur, rhoray, trup, trdown );

	for (ib=0; ib<3; ib++){
		sphalb[ib] = sphalb0[(int)(taur[ib] / TAUSTEP4SPHALB + 0.5)];
		Ttotrayu = ((2 / 3. + muv) + (2 / 3. - muv) * trup[ib])   / (4 / 3. + taur[ib]);
		Ttotrayd = ((2 / 3. + mus) + (2 / 3. - mus) * trdown[ib]) / (4 / 3. + taur[ib]);
		tO3 = tO2 = tH2O = 1;
		if (aO3[ib] != 0) tO3 = expf(-m * UO3 * aO3[ib]);
		if (bH2O[ib] != 0) tH2O = expf(-expf(aH2O[ib] + bH2O[ib] * logf(m * UH2O)));
		TtotraytH2O[ib] = Ttotrayu * Ttotrayd * tH2O;
		tOG[ib] = tO3 * tO2;
	}
	return 0;
}


/// 针对特定波段进行快速大气校正计算
/// @param refl 该波段卫星接收到的大气顶层反射率[0~1.0]
/// @param TtotraytH2O 该波段Rayleigh和H2O综合透过率
/// @param tOG 该波段O2和O3综合透过率
/// @param rhoray 该波段Rayleigh散射反射率
/// @param sphalb 该波段球反照率
/// @return 返回大气校正后的反射率
float correctedrefl(float refl, float TtotraytH2O, float tOG, float rhoray, float sphalb)
{
	float corr_refl;
	corr_refl = (refl / tOG - rhoray) / TtotraytH2O;
	corr_refl /= (1 + corr_refl * sphalb);
	return corr_refl;
}

void linearScale255( float* arr , int asize , short* resarr )
{
	for(int it = 0 ; it<asize ; ++ it )
	{
		resarr[it] = arr[it] * 231.9f ;
	}
}

void jacScale255( short* arr , int asize , short* resarr )
{
	for(int it = 0 ; it<asize ; ++ it )
	{
		if( arr[it] < 30 )
		{
			resarr[it] = arr[it] * 3.667f ;
		}else if( arr[it] < 60 )
		{
			resarr[it] = arr[it] * 1.667f + 60  ;
		}else if( arr[it] < 120 )
		{
			resarr[it] = arr[it] * 0.833f + 110 ;
		}else if( arr[it] < 190 )
		{
			resarr[it] = arr[it] * 0.4286f + 158.57 ;
		}else
		{
			resarr[it] = arr[it] * 0.2308f + 196.15 ;
		}
	}
}


void fy3dgeofile(string b250file, string& geo1kfile, string& geoqkfile)
{
	int len = b250file.length();
	if (len >= 45)
	{
		geo1kfile = wft_replaceString(b250file, "0250M_MS", "GEO1K_MS");
		geoqkfile = wft_replaceString(b250file, "0250M_MS", "GEOQK_MS");
	}
}


int main(int argc, char* argv[])
{
	//
	cout<<"fy3truecolor mod/fy3b/fy3d 0250M.HDF dem.tif outputrgb.tif"<<endl ;

	if( argc != 5 )
	{
		cout<<"out"<<endl ;
		return 11;
	}

	GDALAllRegister() ;

	string type = argv[1] ;
	string outfile = argv[4];

	int xsize = 0;
	int ysize = 0;

	if (type == "fy3d")
	{
		//FY3D_MERSI_GBAL_L1_20180308_0615_0250M_MS.HDF
		//FY3D_MERSI_GBAL_L1_20180308_0615_GEO1K_MS.HDF
		//FY3D_MERSI_GBAL_L1_20180308_0615_GEOQK_MS.HDF
		string b250file = argv[2];
		string demfile = argv[3];
		string geo1kfile, geoqkfile;
		fy3dgeofile(b250file , geo1kfile , geoqkfile );

		if (wft_test_file_exists(geo1kfile) == false)
		{
			cout << "Can not find file " << geo1kfile << endl;
			return 1;
		}
		
		if (wft_test_file_exists(geoqkfile) == false)
		{
			cout << "Can not find file " << geoqkfile << endl;
			return 1;
		}
		
		string band1path = "HDF5:\"" + b250file + "\"://Data/EV_250_RefSB_b1";
		GDALDataset* ds = (GDALDataset*)GDALOpen(band1path.c_str(), GA_ReadOnly);
		xsize = ds->GetRasterXSize();
		ysize = ds->GetRasterYSize();
		GDALClose(ds);

		string p1path = "HDF5:\"" + geo1kfile + "\"://Geolocation/SensorAzimuth";
		string p0path = "HDF5:\"" + geo1kfile + "\"://Geolocation/SolarAzimuth";
		string t1path = "HDF5:\"" + geo1kfile + "\"://Geolocation/SensorZenith";
		string t0path = "HDF5:\"" + geo1kfile + "\"://Geolocation/SolarZenith";

		string t0file = wft_changetail(outfile, ".tif", "_t0.tif");
		string t1file = wft_changetail(outfile, ".tif", "_t1.tif");
		string p0file = wft_changetail(outfile, ".tif", "_p0.tif");
		string p1file = wft_changetail(outfile, ".tif", "_p1.tif");

		string arr0[4] = { t0path  , t1path , p0path , p1path };
		string arr1[4] = { t0file,t1file,p0file,p1file };

		for (int i = 0; i < 4; ++i)
		{
			string cmd1 = gdaltranslate + " -outsize " + wft_int2str(xsize) + " " + wft_int2str(ysize) +
				" -r bilinear " + arr0[i] + " " + arr1[i];
			cout << cmd1 << endl;
			system(cmd1.c_str());
		}
		


		
		return 0;

	}



	// string dir = "E:/coding/MODISTRUECOLOR/fy3d/" ;
	string bandfile = argv[2] ;

	string demfile = argv[3] ;

	string lonfile = argv[8] ;
	string latfile = argv[9] ;

	string t0file = argv[4] ;
	string t1file = argv[6] ;
	string p0file = argv[5] ;
	string p1file = argv[7] ;

	int  asize , demxsize , demysize ;

	float* band1 = 0 ;
	float* band2 = 0 ;
	float* band3 = 0 ;
	float* lonarr = 0 ;
	float* latarr = 0 ;
	float* t0arr = 0 ;
	float* t1arr = 0 ;
	float* p0arr = 0 ;
	float* p1arr = 0 ;
	short* demarr = 0 ;

	float* oband1 = 0 ;
	float* oband2 = 0 ;
	float* oband3 = 0 ;

	float* aband1 = 0 ;
	float* aband2 = 0 ;
	float* aband3 = 0 ;

	///dem
	{
		GDALDataset* ds = (GDALDataset*) GDALOpen( demfile.c_str() , GA_ReadOnly) ;
		demxsize = ds->GetRasterXSize() ;
		demysize = ds->GetRasterYSize() ;
		demarr = new short[demxsize*demysize] ;
		ds->GetRasterBand(1)->RasterIO(GF_Read,0,0,demxsize,demysize,demarr
										,demxsize,demysize,GDT_Int16,0,0,0) ;
		GDALClose(ds) ;
	}

	{
		GDALDataset* ds = (GDALDataset*) GDALOpen( bandfile.c_str() , GA_ReadOnly) ;
		xsize = ds->GetRasterXSize() ;
		ysize = ds->GetRasterYSize() ;

		asize = xsize * ysize ;
		band1 = new float[asize] ;
		band2 = new float[asize] ;
		band3 = new float[asize] ;

		oband1 = new float[asize] ;
		oband2 = new float[asize] ;
		oband3 = new float[asize] ;

		aband1 = new float[asize] ;
		aband2 = new float[asize] ;
		aband3 = new float[asize] ;

		ds->GetRasterBand(1)->RasterIO(GF_Read,0,0,xsize,ysize,band1,xsize,ysize,GDT_Float32,0,0,0) ;
		ds->GetRasterBand(2)->RasterIO(GF_Read,0,0,xsize,ysize,band2,xsize,ysize,GDT_Float32,0,0,0) ;
		ds->GetRasterBand(3)->RasterIO(GF_Read,0,0,xsize,ysize,band3,xsize,ysize,GDT_Float32,0,0,0) ;

		GDALClose(ds) ;
	}

	{
		GDALDataset* ds = (GDALDataset*) GDALOpen( lonfile.c_str() , GA_ReadOnly) ;
		lonarr = new float[asize] ;
		ds->GetRasterBand(1)->RasterIO(GF_Read,0,0,xsize,ysize,lonarr,xsize,ysize,GDT_Float32,0,0,0) ;
		GDALClose(ds) ;
	}

	{
		GDALDataset* ds = (GDALDataset*) GDALOpen( latfile.c_str() , GA_ReadOnly) ;
		latarr = new float[asize] ;
		ds->GetRasterBand(1)->RasterIO(GF_Read,0,0,xsize,ysize,latarr,xsize,ysize,GDT_Float32,0,0,0) ;
		GDALClose(ds) ;
	}

	{
		GDALDataset* ds = (GDALDataset*) GDALOpen( t0file.c_str() , GA_ReadOnly) ;
		t0arr = new float[asize] ;
		ds->GetRasterBand(1)->RasterIO(GF_Read,0,0,xsize,ysize,t0arr,xsize,ysize,GDT_Float32,0,0,0) ;
		GDALClose(ds) ;
	}

	{
		GDALDataset* ds = (GDALDataset*) GDALOpen( t1file.c_str() , GA_ReadOnly) ;
		t1arr = new float[asize] ;
		ds->GetRasterBand(1)->RasterIO(GF_Read,0,0,xsize,ysize,t1arr,xsize,ysize,GDT_Float32,0,0,0) ;
		GDALClose(ds) ;
	}

	{
		GDALDataset* ds = (GDALDataset*) GDALOpen( p0file.c_str() , GA_ReadOnly) ;
		p0arr = new float[asize] ;
		ds->GetRasterBand(1)->RasterIO(GF_Read,0,0,xsize,ysize,p0arr,xsize,ysize,GDT_Float32,0,0,0) ;
		GDALClose(ds) ;
	}

	{
		GDALDataset* ds = (GDALDataset*) GDALOpen( p1file.c_str() , GA_ReadOnly) ;
		p1arr = new float[asize] ;
		ds->GetRasterBand(1)->RasterIO(GF_Read,0,0,xsize,ysize,p1arr,xsize,ysize,GDT_Float32,0,0,0) ;
		GDALClose(ds) ;
	}

	double a2arr[] = { 0 , 0 , 0 } ;
	double a1arr[] = { 0.02542 , 0.02682 ,  0.02783 } ;
	double a0arr[] = { -3.253,-3.609,-6.913 } ;

	if( type == "mod" )
	{
		a2arr[0] = a2arr[1] = a2arr[2] = 0 ;
		a1arr[0] = a1arr[1] = a1arr[2] = 1 ;
		a0arr[0] = a0arr[1] = a0arr[2] = 0 ;
	}else if( type == "fy3b" )
	{
		a2arr[0] = a2arr[1] = a2arr[2] = 0 ;
		a1arr[0] = a1arr[1] = a1arr[2] = 0.0001 ;
		a0arr[0] = a0arr[1] = a0arr[2] = 0 ;
	}else if( type=="fy3d" )
	{
		a2arr[0] = a2arr[1] = a2arr[2] = 0 ;
		a1arr[0] = 0.02542 * 0.01 ;
		a1arr[1] = 0.02682 * 0.01 ;
		a1arr[2] = 0.02783 * 0.01 ;

		a0arr[0] = -3.253 * 0.01 ;
		a0arr[1] = -3.609 * 0.01 ;
		a0arr[2] = -6.913 * 0.01 ;
	}


	float* bandPtrArr[] = { band1 , band2, band3} ;
	float* obandPtrArr[] = { oband1 , oband2, oband3} ;
	float* abandPtrArr[] = { aband1 , aband2, aband3} ;
	cout<<"Processing..."<<endl ;
	int percentOld = 0 ;
	for(int it = 0 ; it < asize ; ++ it )
	{
		short height = getDemByLonLat( demarr , demxsize , demysize , lonarr[it] , latarr[it] ) ;
		float t0 = t0arr[it] * 0.01 ;
		float t1 = t1arr[it] * 0.01 ;
		float ph0 = p0arr[it] * 0.01 ;
		float ph1 = p1arr[it] * 0.01 ;

		float phi = ph0 - ph1 ;
		float mus = cos( t0*DEG2RAD) ;
		float muv = cos( t1*DEG2RAD) ;

		float sphalb[3]={0,0,0} ;
		float rhoray[3] = {0,0,0} ;
		float TtotraytH2O[3] = {0,0,0} ;
		float tOG[3] = {0,0,0} ;
		float F0[3] = {1877,1862,1606} ;
		int atmok = getatmvariables(
				 mus, 
				 muv, 
				 phi, 
				 height, 
				sphalb, 
				rhoray, 
				TtotraytH2O, 
				tOG
				) ;
		if( atmok == 0 )
		{//ok
			for(int iband = 0 ; iband<3 ; ++ iband )
			{
				float* bandPtr = bandPtrArr[iband] ;
				float* obandPtr = obandPtrArr[iband] ;
				float* abandPtr = abandPtrArr[iband] ;

				double bval = bandPtr[it]*1.0 ;
				double refl = bval*bval*a2arr[iband] + bval*a1arr[iband] + a0arr[iband] ;

				abandPtr[it] = refl ;
				double rtoa = refl/mus ;
				
				obandPtr[it] = rtoa ;

				double acVal = correctedrefl( rtoa , 
					 TtotraytH2O[iband] ,
					 tOG[iband] ,  
					 rhoray[iband] ,  
					 sphalb[iband] ) ;
				bandPtr[it] = acVal ;
			}
		}else
		{//bad atm var
			band1[it] = -0.01 ;
			band2[it] = -0.01 ;
			band3[it] = -0.01 ;
		}	
		if( it % 1000 == 0 ){
			int percent = it*1.f/asize * 100 +0.5  ;
			if( percent%5==0 && percent != percentOld )
			{
				percentOld = percent ;
				cout<<percent<<".." ;
			}
		}
	}
	cout<<endl ;

	{
		string outb3file = wft_base_name(bandfile) + "-ac.tif" ;
		GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset* dsout = driver->Create(outb3file.c_str() , xsize,ysize,3,GDT_Float32,0) ;
		dsout->GetRasterBand(1)->RasterIO(GF_Write,0,0,xsize,ysize,band1,xsize,ysize,GDT_Float32,0,0,0) ;
		dsout->GetRasterBand(2)->RasterIO(GF_Write,0,0,xsize,ysize,band2,xsize,ysize,GDT_Float32,0,0,0) ;
		dsout->GetRasterBand(3)->RasterIO(GF_Write,0,0,xsize,ysize,band3,xsize,ysize,GDT_Float32,0,0,0) ;
		GDALClose(dsout) ;	

		string cmd1 = "gdal_translate " + outb3file + " " + outb3file + ".png" ;
		system(cmd1.c_str()) ;
	}

	{
		string outb3file = wft_base_name(bandfile) + "-acbyte.tif" ;
		GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset* dsout = driver->Create(outb3file.c_str() , xsize,ysize,3,GDT_Byte,0) ;
		short* resarr = new short[asize] ;
		linearScale255(band1 , asize , resarr) ;
		dsout->GetRasterBand(3)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		linearScale255(band2 , asize , resarr) ;
		dsout->GetRasterBand(2)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		linearScale255(band3 , asize , resarr) ;
		dsout->GetRasterBand(1)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		GDALClose(dsout) ;

		string cmd1 = "gdal_translate " + outb3file + " " + outb3file + ".png" ;
		system(cmd1.c_str()) ;

		delete [] resarr ;
	}

	{
		string outb3file = wft_base_name(bandfile) + "-toa.tif" ;
		GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset* dsout = driver->Create(outb3file.c_str() , xsize,ysize,3,GDT_Byte,0) ;
		short* resarr = new short[asize] ;
		linearScale255(oband1 , asize , resarr) ;
		dsout->GetRasterBand(3)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		linearScale255(oband2 , asize , resarr) ;
		dsout->GetRasterBand(2)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		linearScale255(oband3 , asize , resarr) ;
		dsout->GetRasterBand(1)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		GDALClose(dsout) ;	

		string cmd1 = "gdal_translate " + outb3file + " " + outb3file + ".png" ;
		system(cmd1.c_str()) ;

		delete [] resarr ;
	}

	{
		string outb3file = wft_base_name(bandfile) + "-refl.tif" ;
		GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset* dsout = driver->Create(outb3file.c_str() , xsize,ysize,3,GDT_Byte,0) ;
		short* resarr = new short[asize] ;
		linearScale255(aband1 , asize , resarr) ;
		dsout->GetRasterBand(3)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		linearScale255(aband2 , asize , resarr) ;
		dsout->GetRasterBand(2)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		linearScale255(aband3 , asize , resarr) ;
		dsout->GetRasterBand(1)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		GDALClose(dsout) ;	

		string cmd1 = "gdal_translate " + outb3file + " " + outb3file + ".png" ;
		system(cmd1.c_str()) ;

		delete [] resarr ;
	}

	{
		string outb3file = wft_base_name(bandfile) + "-acjac.tif" ;
		GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset* dsout = driver->Create(outb3file.c_str() , xsize,ysize,3,GDT_Byte,0) ;
		short* resarr = new short[asize] ;
		linearScale255( band1 , asize , resarr) ;
		jacScale255( resarr , asize , resarr) ;
		dsout->GetRasterBand(3)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		linearScale255( band2 , asize , resarr) ;
		jacScale255( resarr , asize , resarr) ;
		dsout->GetRasterBand(2)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		linearScale255( band3 , asize , resarr) ;
		jacScale255( resarr , asize , resarr) ;
		dsout->GetRasterBand(1)->RasterIO(GF_Write,0,0,xsize,ysize,resarr,xsize,ysize,GDT_Int16,0,0,0) ;
		GDALClose(dsout) ;	

		string cmd1 = "gdal_translate " + outb3file + " " + outb3file + ".png" ;
		system(cmd1.c_str()) ;

		delete [] resarr ;
	}

	delete[] band1;
	delete[] band2 ;
	delete[] band3 ;

	delete[] aband1;
	delete[] aband2 ;
	delete[] aband3 ;

	delete[] oband1;
	delete[] oband2 ;
	delete[] oband3 ;

	delete[] lonarr ; delete [] latarr ;
	delete[] t0arr ; delete[] t1arr ; delete[] p0arr; delete[] p1arr ;
	
	cout<<"done."<<endl ;
	return 0;
}

