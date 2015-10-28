function sinh(x) { return 0.5*(Math.exp(x)-Math.exp(-x)) }
function cosh(x) { return 0.5*(Math.exp(x)+Math.exp(-x)) }
function tanh(x) { return sinh(x)/cosh(x)}
function atanh(x) { return Math.log((1+x)/(1-x)) }
function arctanh(x) { return 0.5*Math.log((1+x)/(1-x)) } 



// 地球の赤道半径
a=6378137;

rf=298.257222101;https://github.com/Kenya39/coordinate/edit/test/lib/coordinate.js#
m0=0.9999;
s2r=Math.PI/648000;
n=0.5/(rf-0.5);
n15=1.5*n ;
anh=0.5*a/(1+n) ;
nsq=n*n ;
e2n=2*Math.sqrt(n)/(1+n) ;
ra=2*anh*m0*(1+nsq/4+nsq*nsq/64)
jt=5 ; jt2=2*jt ; ep=1.0 ; e=[] ; s=[0.0] ; t=[] ;
alp=[]
beta=[] ; dlt=[]

 for(k=1; k<=jt; k++) { ep*=e[k]=n15/k-n ; e[k+jt]=n15/(k+jt)-n }

// 展開パラメータの事前入力
alp[1]=(1/2+(-2/3+(5/16+(41/180-127/288*n)*n)*n)*n)*n
alp[2]=(13/48+(-3/5+(557/1440+281/630*n)*n)*n)*nsq 
alp[3]=(61/240+(-103/140+15061/26880*n)*n)*n*nsq
alp[4]=(49561/161280-179/168*n)*nsq*nsq
alp[5]=34729/80640*n*nsq*nsq  

// 展開パラメータの事前入力
 beta[1]=(1/2+(-2/3+(37/96+(-1/360-81/512*n)*n)*n)*n)*n
 beta[2]=(1/48+(1/15+(-437/1440+46/105*n)*n)*n)*nsq
 beta[3]=(17/480+(-37/840-209/4480*n)*n)*n*nsq
 beta[4]=(4397/161280-11/504*n)*nsq*nsq
 beta[5]=4583/161280*n*nsq*nsq
 dlt[1]=(2+(-2/3+(-2+(116/45+(26/45-2854/675*n)*n)*n)*n)*n)*n
 dlt[2]=(7/3+(-8/5+(-227/45+(2704/315+2323/945*n)*n)*n)*n)*nsq
 dlt[3]=(56/15+(-136/35+(-1262/105+73814/2835*n)*n)*n)*n*nsq
 dlt[4]=(4279/630+(-332/35-399572/14175*n)*n)*nsq*nsq
 dlt[5]=(4174/315-144838/6237*n)*n*nsq*nsq
 dlt[6]=601676/22275*nsq*nsq*nsq

// 平面直角座標の座標系原点の緯度を度単位で、経度を分単位で格納
 phi0=[0,33,33,36,33,36,36,36,36,36,40,44,44,44,26,26,26,26,20,26]
 lmbd0=[0,7770,7860,7930,8010,8060,8160,8230,8310,8390,8450,8415,8535,8655,8520,7650,7440,7860,8160,9240]

// 該当緯度の 2 倍角の入力により赤道からの子午線弧長を求める関数
function Merid(phi2) { 
   dc=2.0*Math.cos(phi2) ; s[1]=Math.sin(phi2)
   for(i=1; i<=jt2; i++) { s[i+1]=dc*s[i]-s[i-1] ; t[i]=(1.0/i-4.0*i)*s[i] }
   sum=0.0 ; c1=ep ; j=jt
   while(j) {
     c2=phi2 ; c3=2.0 ; l=j ; m=0
     while(l) { c2+=(c3/=e[l--])*t[++m]+(c3*=e[2*j-l])*t[++m] }
     sum+=c1*c1*c2 ; c1/=e[j--]
   }
   return anh*(sum+phi2)
}  

// 緯度Lat,経度Lonの設定
var lon = 0.0; var lat = 0.0;

// ズームレベル
var zl = 17;

// ピクセル座標変換時の定数L
// L = 180/π * asin(tanh(π))
var L = (180/Math.PI)*(Math.asin(tanh(Math.PI)));



function lon2tile(lon,zoom) {

  return (Math.floor((lon+180)/360*Math.pow(2,zoom))); }

//function lat2tile(lat,zoom)  { return (Math.floor((1-Math.log(Math.tan(lat*Math.PI/180) + 1/Math.cos(lat*Math.PI/180))/Math.PI)/2 *Math.pow(2,zoom))); }
function lat2tile(lat,zoom)  {
  var a = lat*Math.PI/180;
  var b = Math.abs(Math.tan(a) + 1/Math.cos(a));
  var tile = (1-Math.log(b)/Math.PI) * Math.pow(2,zoom-1)
  var tileFl = Math.floor(tile);
  return (tileFl);
}
 
function tile2lon(x,zoom) {
  return (x/Math.pow(2,zoom)*360-180);
}
function tile2lat(y,zoom) {
  var n=Math.PI-2*Math.PI*y/Math.pow(2,zoom);
  return (180/Math.PI*Math.atan(0.5*(Math.exp(n)-Math.exp(-n))));
}

function tile2extent(tileZXY){
  var minLon = tile2lon( tileZXY[1],tileZXY[0] );
  var minLat = tile2lat( tileZXY[2]+1, tileZXY[0] );
  var maxLon = tile2lon( tileZXY[1]+1,tileZXY[0] );
  var maxLat = tile2lat( tileZXY[2],tileZXY[0] );
  
  var extent = [minLon, minLat, maxLon, maxLat];
  
  return extent;
}


// 経度、緯度からピクセル座標X、Yを変換する
function LonLat2PxXPxY(zl,Lon,Lat){
    // 緯度経度からピクセル座標を求める
    // 参考：http://www.trail-note.net/tech/coordinate/
    var pxX = Math.pow(2,zl+7)*(lon/180+1)
//    var pxY = (Math.pow(2,zl+7)/Math.PI)*(-1*atanh(Math.sin((Math.PI/180)*lon))+atanh(Math.sin((Math.PI/180)*L)))
//  国土地理院用？
    var pxY = (Math.pow(2,zl+6)/Math.PI)*(-1*atanh(Math.sin((Math.PI/180)*lat))+atanh(Math.sin((Math.PI/180)*L)))

    return [pxX,pxY];
}

// 座標クラスの生成
function coordinate(kei){

  // 座標系を登録
  this.kei = kei;

  this.getLonLat = function(){

    return this.pointLonLat;

  };


  // 平面直角座標の設定
  // 平面直角座標が設定された時点で緯度経度を計算して登録する
  this.setXY = function(x,y){
  
    this.X = x; this.Y = y;
    this.pointXY = [x,y];

    // 実際の計算実行部分
    xip=xi=(y+m0*Merid(2*phi0[kei]*3600*s2r))/ra ; etap=eta=x/ra ; sgmp=1 ; taup=0
      for(j=beta.length; --j; ) {
        besin=beta[j]*Math.sin(2*j*xi) ; becos=beta[j]*Math.cos(2*j*xi)
        xip-=besin*cosh(2*j*eta) ; etap-=becos*sinh(2*j*eta)
        sgmp-=2*j*becos*cosh(2*j*eta) ; taup+=2*j*besin*sinh(2*j*eta)
      }
    sxip=Math.sin(xip) ; cxip=Math.cos(xip) ; shetap=sinh(etap) ; chetap=cosh(etap)
    phi=chi=Math.asin(sxip/chetap)
    for(j=dlt.length; --j; ) { phi+=dlt[j]*Math.sin(2*j*chi) }
    nphi=(1-n)/(1+n)*Math.tan(phi)

    lmbd=lmbd0[kei]*60+Math.atan2(shetap, cxip)/s2r
    gmm=Math.atan2(taup*cxip*chetap+sgmp*sxip*shetap, sgmp*cxip*chetap-taup*sxip*shetap)
    m=ra/a*Math.sqrt((cxip*cxip+shetap*shetap)/(sgmp*sgmp+taup*taup)*(1+nphi*nphi))
 
    // ラジアン → 度分秒変換
    ido=Math.floor((phi/s2r/3600)+0.0000000001)
    ifun=Math.floor(((phi/s2r-ido*3600)/60)+0.0000000001)
    ibyou=Math.floor((phi/s2r-ido*3600-ifun*60)*100000+0.5)/100000
 
    keido=Math.floor((lmbd/3600)+0.0000000001)
    keifun=Math.floor(((lmbd-keido*3600)/60)+0.0000000001)
    keibyou=Math.floor((lmbd-keido*3600-keifun*60)*100000+0.5)/100000
 
    sgn=(gmm<0)
    gdo=Math.floor(gmm/s2r/3600)+sgn
    gfun=Math.floor((gmm/s2r-gdo*3600)/60)+sgn
    gbyou=gmm/s2r-gdo*3600-gfun*60


    // 緯度Lat,経度Lonの設定
    lon = Math.floor((lmbd/3600)*10000000000+0.01)/10000000000 ;
    lat = Math.floor((phi/s2r/3600)*10000000000+0.01)/10000000000 ;
    this.Lon = lon; this.Lat = lat;
    this.pointLonLat = [lon,lat];


// 結果表示

//    document.write("<h2>座標系番号： " + kei + "  Ｘ座標： " + x + "  Ｙ座標： " + y + "<br/><br/>") 

//    document.write("phi＝" + phi + "″<br/>")
//    document.write("phi/s2r＝" + phi/s2r + "″<br/>")
//    document.write("φ＝" + ido + "°" + ifun + "′" + ibyou + "″，λ＝" + keido + "° " + keifun + "′" + keibyou + "″<br/>")
//    document.write("lmbd＝" + lmbd + "″<br/>")
//    document.write("γ＝" + (sgn?"－":"＋") + Math.abs(gdo) + "°" + Math.abs(gfun) + "′ " + Math.abs(gbyou) + "″，m＝" + m + "<br/></h2>")  

  
  };
  // 対象のグーグルタイルを取得
  // 1 zl  ズームレベル
  // 2 lon 経度
  // 2 lat 緯度
  this.getGTileZXY = function(zl){

    var pointPix = LonLat2PxXPxY(zl, lon, lat);
    
    var gtileX = Math.floor((pointPix[0] / 256)+0.0000000001);
    var gtileY = Math.floor((pointPix[1] / 256)+0.0000000001);

    var gtileZXY = [zl,gtileX,gtileY]

// sample
//    document.write("http://cyberjapandata.gsi.go.jp/xyz/std/" + zl + "/" + gtileX + "/" + gtileY +".png<br/>")

    return gtileZXY;
    
    

  };
  
  // 対象のSlippy Map tilenameを取得
  // 1 zl  ズームレベル
  // 2 lon 経度
  // 2 lat 緯度
  this.getSTileZXY = function(zl){

   
    var stileX = lon2tile(lon,zl);
    var stileY = lat2tile(lat,zl);

    var stileZXY = [zl,stileX,stileY]

// sample
//   document.write("http://a.tile.openstreetmap.org/" + zl + "/" + stileX + "/" + stileY +".png<br/>")

    return stileZXY;    

  };
  // 経度緯度設定して、平面直角座標を計算する
  this.setLonLat = function( Lon, Lat ){

    // 与件入力
//  num=eval(prompt("座標系番号を入力してください。"))
//  phi=eval(prompt("緯度を ddmmss.ssss 形式で入力してください。"))
    latdeg=Math.floor(Lat+0.00000000001) ;
    latminOrg=(Lat-latdeg)*60;
    latmin=Math.floor(latminOrg+0.000000001);
    latsec=(latminOrg-latmin)*60
    phi = latdeg*10000 + latmin*100 + latsec;


//  lmbd=eval(prompt("経度を dddmmss.ssss 形式で入力してください。"))

    phideg=Math.floor(phi/10000+0.00000000001) ;
    phimin=Math.floor((phi-phideg*10000)/100+0.000000001)
    phirad=(phideg*3600+phimin*60+phi-phideg*10000-phimin*100)*s2r
//    alert(phirad);
//    phirad1 = Lat * Math.PI / 180 ;
//    alert(phirad1);

    londeg=Math.floor(Lon+0.00000000001) ;
    lonminOrg=(Lon-londeg)*60
    lonmin=Math.floor(lonminOrg+0.000000001)
    lonsec=(lonminOrg-lonmin)*60
    lmbd = londeg*10000 + lonmin*100 + lonsec;
//    alert(lmbd);

    num = kei;

    lmbddeg=Math.floor(lmbd/10000+0.00000000001) ;
    lmbdmin=Math.floor((lmbd-lmbddeg*10000)/100+0.000000001) 
    lmbdsec=lmbddeg*3600+lmbdmin*60+lmbd-lmbddeg*10000-lmbdmin*100 
 
    // 実際の計算実行部分 
    sphi=Math.sin(phirad) ; nphi=(1-n)/(1+n)*Math.tan(phirad) 
    dlmbd=(lmbdsec-lmbd0[num]*60)*s2r 
    sdlmbd=Math.sin(dlmbd) ; cdlmbd=Math.cos(dlmbd) 
    tchi=sinh(arctanh(sphi)-e2n*arctanh(e2n*sphi)) ; cchi=Math.sqrt(1+tchi*tchi)
    
    xi=xip=Math.atan2(tchi, cdlmbd) ;

    eta=etap=arctanh(sdlmbd/cchi) ; sgm=1 ; tau=0 
    for(j=alp.length; --j; ) { 
        alsin=alp[j]*Math.sin(2*j*xip) ; alcos=alp[j]*Math.cos(2*j*xip) 
        xi+=alsin*cosh(2*j*etap) ; eta+=alcos*sinh(2*j*etap) 
        sgm+=2*j*alcos*cosh(2*j*etap) ; tau+=2*j*alsin*sinh(2*j*etap) 
    } 
 
    x=ra*eta;

    y=ra*xi-m0*Merid(2*phi0[num]*3600*s2r) ;

    
    this.X = x;
    this.Y = y;
    gmm=Math.atan2(tau*cchi*cdlmbd+sgm*tchi*sdlmbd, sgm*cchi*cdlmbd-tau*tchi*sdlmbd) 
    m=ra/a*Math.sqrt((sgm*sgm+tau*tau)/(tchi*tchi+cdlmbd*cdlmbd)*(1+nphi*nphi)) 

    // ラジアン → 度分秒変換 
    sgn=(gmm<0) 
    gdo=Math.floor(gmm/s2r/3600)+sgn 
    gfun=Math.floor((gmm/s2r-gdo*3600)/60)+sgn 
    gbyou=gmm/s2r-gdo*3600-gfun*60 

    // 結果表示
//    document.write("<h2>座標系番号： " + kei + "  緯度： " + phi + "  経度： " + lmbd + "<br/><br/>")
//    document.write("Ｘ＝" + x + "，Ｙ＝" + y + "<br/>")
//    document.write("γ＝" + (sgn?"－":"＋") + Math.abs(gdo) + "°" + Math.abs(gfun) + "′ " + Math.abs(gbyou) + "″，m＝" + m + "<br/></h2>")  
  
  };
}

