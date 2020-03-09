     var forEachPoint = function(func) {
        return function(input, opt_output, opt_dimension) {
          var len = input.length;
          var dimension = opt_dimension ? opt_dimension : 2;
          var output;
          if (opt_output) {
            output = opt_output;
          } else {
            if (dimension !== 2) {
              output = input.slice();
            } else {
              output = new Array(len);
            }
          }
          for (var offset = 0; offset < len; offset += dimension) {
            func(input, output, offset)
          }
          return output;
        };
      };
      
      var sphericalMercator = {}
      
      var RADIUS = 6378137;
      var MAX_LATITUDE = 85.0511287798;
      var RAD_PER_DEG = Math.PI / 180;
      
      sphericalMercator.forward = forEachPoint(function(input, output, offset) {
        var lat = Math.max(Math.min(MAX_LATITUDE, input[offset + 1]), -MAX_LATITUDE);
        var sin = Math.sin(lat * RAD_PER_DEG);
      
        output[offset] = RADIUS * input[offset] * RAD_PER_DEG;
        output[offset + 1] = RADIUS * Math.log((1 + sin) / (1 - sin)) / 2;
      });
      
      sphericalMercator.inverse = forEachPoint(function(input, output, offset) {
        output[offset] = input[offset] / RADIUS / RAD_PER_DEG;
        output[offset + 1] = (2 * Math.atan(Math.exp(input[offset + 1] / RADIUS)) - (Math.PI / 2)) / RAD_PER_DEG;
      });
      
      
      var baiduMercator = {}
      
      var MCBAND = [12890594.86, 8362377.87,
          5591021, 3481989.83, 1678043.12, 0];
      
      var LLBAND = [75, 60, 45, 30, 15, 0];
      
      var MC2LL = [
          [1.410526172116255e-8, 0.00000898305509648872, -1.9939833816331,
              200.9824383106796, -187.2403703815547, 91.6087516669843,
              -23.38765649603339, 2.57121317296198, -0.03801003308653,
              17337981.2],
          [-7.435856389565537e-9, 0.000008983055097726239,
              -0.78625201886289, 96.32687599759846, -1.85204757529826,
              -59.36935905485877, 47.40033549296737, -16.50741931063887,
              2.28786674699375, 10260144.86],
          [-3.030883460898826e-8, 0.00000898305509983578, 0.30071316287616,
              59.74293618442277, 7.357984074871, -25.38371002664745,
              13.45380521110908, -3.29883767235584, 0.32710905363475,
              6856817.37],
          [-1.981981304930552e-8, 0.000008983055099779535, 0.03278182852591,
              40.31678527705744, 0.65659298677277, -4.44255534477492,
              0.85341911805263, 0.12923347998204, -0.04625736007561,
              4482777.06],
          [3.09191371068437e-9, 0.000008983055096812155, 0.00006995724062,
              23.10934304144901, -0.00023663490511, -0.6321817810242,
              -0.00663494467273, 0.03430082397953, -0.00466043876332,
              2555164.4],
          [2.890871144776878e-9, 0.000008983055095805407, -3.068298e-8,
              7.47137025468032, -0.00000353937994, -0.02145144861037,
              -0.00001234426596, 0.00010322952773, -0.00000323890364,
              826088.5]];
      
      var LL2MC = [
          [-0.0015702102444, 111320.7020616939, 1704480524535203,
              -10338987376042340, 26112667856603880,
              -35149669176653700, 26595700718403920,
              -10725012454188240, 1800819912950474, 82.5],
          [0.0008277824516172526, 111320.7020463578, 647795574.6671607,
              -4082003173.641316, 10774905663.51142, -15171875531.51559,
              12053065338.62167, -5124939663.577472, 913311935.9512032,
              67.5],
          [0.00337398766765, 111320.7020202162, 4481351.045890365,
              -23393751.19931662, 79682215.47186455, -115964993.2797253,
              97236711.15602145, -43661946.33752821, 8477230.501135234,
              52.5],
          [0.00220636496208, 111320.7020209128, 51751.86112841131,
              3796837.749470245, 992013.7397791013, -1221952.21711287,
              1340652.697009075, -620943.6990984312, 144416.9293806241,
              37.5],
          [-0.0003441963504368392, 111320.7020576856, 278.2353980772752,
              2485758.690035394, 6070.750963243378, 54821.18345352118,
              9540.606633304236, -2710.55326746645, 1405.483844121726,
              22.5],
          [-0.0003218135878613132, 111320.7020701615, 0.00369383431289,
              823725.6402795718, 0.46104986909093, 2351.343141331292,
              1.58060784298199, 8.77738589078284, 0.37238884252424, 7.45]];
      
      
      function getRange(v, min, max) {
        v = Math.max(v, min);
        v = Math.min(v, max);
      
        return v;
      }
      
      function getLoop(v, min, max) {
        var d = max - min;
        while (v > max) {
          v -= d;
        }
        while (v < min) {
          v += d;
        }
      
        return v;
      }
      
      function convertor(input, output, offset, table) {
        var px = input[offset];
        var py = input[offset + 1];
        var x = table[0] + table[1] * Math.abs(px);
        var d = Math.abs(py) / table[9];
        var y = table[2]
            + table[3]
            * d
            + table[4]
            * d
            * d
            + table[5]
            * d
            * d
            * d
            + table[6]
            * d
            * d
            * d
            * d
            + table[7]
            * d
            * d
            * d
            * d
            * d
            + table[8]
            * d
            * d
            * d
            * d
            * d
            * d;
      
        output[offset] = x * (px < 0 ? -1 : 1);
        output[offset + 1] = y * (py < 0 ? -1 : 1);
      }
      
      baiduMercator.forward = forEachPoint(function(input, output, offset) {
        var lng = getLoop(input[offset], -180, 180);
        var lat = getRange(input[offset + 1], -74, 74);
      
        var table = null;
        var j;
        for (j = 0; j < LLBAND.length; ++j) {
          if (lat >= LLBAND[j]) {
            table = LL2MC[j];
            break;
          }
        }
        if (table === null) {
          for (j = LLBAND.length - 1; j >= 0; --j) {
            if (lat <= -LLBAND[j]) {
              table = LL2MC[j];
              break;
            }
          }
        }
        output[offset] = lng;
        output[offset + 1] = lat;
        convertor(output, output, offset, table);
      });
      
      baiduMercator.inverse = forEachPoint(function(input, output, offset) {
        var y_abs = Math.abs(input[offset + 1]);
      
        var table = null;
        for (var j = 0; j < MCBAND.length; j++) {
          if (y_abs >= MCBAND[j]) {
            table = MC2LL[j];
            break;
          }
        }
      
        convertor(input, output, offset, table);
      });
      
      var gcj02 = {}
      
      var PI = Math.PI;
      var AXIS = 6378245.0;
      var OFFSET = 0.00669342162296594323;  // (a^2 - b^2) / a^2
      
      function delta(wgLon, wgLat) {
        var dLat = transformLat(wgLon - 105.0, wgLat - 35.0);
        var dLon = transformLon(wgLon - 105.0, wgLat - 35.0);
        var radLat = wgLat / 180.0 * PI;
        var magic = Math.sin(radLat);
        magic = 1 - OFFSET * magic * magic;
        var sqrtMagic = Math.sqrt(magic);
        dLat = (dLat * 180.0) / ((AXIS * (1 - OFFSET)) / (magic * sqrtMagic) * PI);
        dLon = (dLon * 180.0) / (AXIS / sqrtMagic * Math.cos(radLat) * PI);
        return [dLon, dLat];
      }
      
      function outOfChina(lon, lat) {
        if (lon < 72.004 || lon > 137.8347) {
          return true;
        }
        if (lat < 0.8293 || lat > 55.8271) {
          return true;
        }
        return false;
      }
      
      function transformLat(x, y) {
        var ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * Math.sqrt(Math.abs(x));
        ret += (20.0 * Math.sin(6.0 * x * PI) + 20.0 * Math.sin(2.0 * x * PI)) * 2.0 / 3.0;
        ret += (20.0 * Math.sin(y * PI) + 40.0 * Math.sin(y / 3.0 * PI)) * 2.0 / 3.0;
        ret += (160.0 * Math.sin(y / 12.0 * PI) + 320 * Math.sin(y * PI / 30.0)) * 2.0 / 3.0;
        return ret;
      }
      
      function transformLon(x, y) {
        var ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * Math.sqrt(Math.abs(x));
        ret += (20.0 * Math.sin(6.0 * x * PI) + 20.0 * Math.sin(2.0 * x * PI)) * 2.0 / 3.0;
        ret += (20.0 * Math.sin(x * PI) + 40.0 * Math.sin(x / 3.0 * PI)) * 2.0 / 3.0;
        ret += (150.0 * Math.sin(x / 12.0 * PI) + 300.0 * Math.sin(x / 30.0 * PI)) * 2.0 / 3.0;
        return ret;
      }
      
      gcj02.toWGS84 = forEachPoint(function(input, output, offset) {
        var lng = input[offset];
        var lat = input[offset + 1];
        if (!outOfChina(lng, lat)) {
          var deltaD = delta(lng, lat);
          lng = lng - deltaD[0];
          lat = lat - deltaD[1];
        }
        output[offset] = lng;
        output[offset + 1] = lat;
      });
      
      gcj02.fromWGS84 = forEachPoint(function(input, output, offset) {
        var lng = input[offset];
        var lat = input[offset + 1];
        if (!outOfChina(lng, lat)) {
          var deltaD = delta(lng, lat);
          lng = lng + deltaD[0];
          lat = lat + deltaD[1];
        }
        output[offset] = lng;
        output[offset + 1] = lat;
      });
      
      var bd09 = {}
      
      var PI = Math.PI;
      var X_PI = PI * 3000 / 180;
      
      function toGCJ02(input, output, offset) {
        var x = input[offset] - 0.0065;
        var y = input[offset + 1] - 0.006;
        var z = Math.sqrt(x * x + y * y) - 0.00002 * Math.sin(y * X_PI);
        var theta = Math.atan2(y, x) - 0.000003 * Math.cos(x * X_PI);
        output[offset] = z * Math.cos(theta);
        output[offset + 1] = z * Math.sin(theta);
        return output;
      }
      
      function fromGCJ02(input, output, offset) {
        var x = input[offset];
        var y = input[offset + 1];
        var z = Math.sqrt(x * x + y * y) + 0.00002 * Math.sin(y * X_PI);
        var theta = Math.atan2(y, x) + 0.000003 * Math.cos(x * X_PI);
        output[offset] = z * Math.cos(theta) + 0.0065;
        output[offset + 1] = z * Math.sin(theta) + 0.006;
        return output;
      }
      
      bd09.toWGS84 = function(input, opt_output, opt_dimension) {
        var output = forEachPoint(toGCJ02)(input, opt_output, opt_dimension);
        return gcj02.toWGS84(output, output, opt_dimension);
      };
      
      bd09.fromWGS84 = function(input, opt_output, opt_dimension) {
        var output = gcj02.fromWGS84(input, opt_output, opt_dimension);
        return forEachPoint(fromGCJ02)(output, output, opt_dimension);
      };
      
      
      var projzh = {}
      
      projzh.smerc2bmerc = function(input, opt_output, opt_dimension) {
        var output = sphericalMercator.inverse(input, opt_output, opt_dimension);
        output = bd09.fromWGS84(output, output, opt_dimension);
        return baiduMercator.forward(output, output, opt_dimension);
      };
      
      projzh.bmerc2smerc = function(input, opt_output, opt_dimension) {
        var output = baiduMercator.inverse(input, opt_output, opt_dimension);
        output = bd09.toWGS84(output, output, opt_dimension);
        return sphericalMercator.forward(output, output, opt_dimension);
      };
      
      projzh.bmerc2ll = function(input, opt_output, opt_dimension) {
        var output = baiduMercator.inverse(input, opt_output, opt_dimension);
        return bd09.toWGS84(output, output, opt_dimension);
      };
      
      projzh.ll2bmerc = function(input, opt_output, opt_dimension) {
        var output = bd09.fromWGS84(input, opt_output, opt_dimension);
        return baiduMercator.forward(output, output, opt_dimension);
      };
      
      projzh.ll2smerc = sphericalMercator.forward;
      projzh.smerc2ll = sphericalMercator.inverse;
      
      
      
      var extent = [72.004, 0.8293, 137.8347, 55.8271];
      
      var baiduMercatorProj = new ol.proj.Projection({
        code: 'baidu',
        extent: ol.extent.applyTransform(extent, projzh.ll2bmerc),
        units: 'm'
      });
      
      ol.proj.addProjection(baiduMercatorProj);
      ol.proj.addCoordinateTransforms('EPSG:4326', baiduMercatorProj, projzh.ll2bmerc, projzh.bmerc2ll);
      ol.proj.addCoordinateTransforms('EPSG:3857', baiduMercatorProj, projzh.smerc2bmerc, projzh.bmerc2smerc);
      
      var bmercResolutions = new Array(19);
      for (var i = 0; i < 19; ++i) {
        bmercResolutions[i] = Math.pow(2, 18 - i);
      }
      var baidu = new ol.layer.Tile({
        source: new ol.source.XYZ({
          projection: 'baidu',
          maxZoom: 18,
          tileUrlFunction: function(tileCoord) {
            var x = tileCoord[1];
            var y = tileCoord[2];
            var z = tileCoord[0];
            return "http://online3.map.bdimg.com/onlinelabel/?qt=tile&x="+x+"&y="+y+"&z="+z+"&styles=pl&udt=20151021&scaler=1&p=1";;
          },
          tileGrid: new ol.tilegrid.TileGrid({
            resolutions: bmercResolutions,
            origin: [0, 0],
            extent: ol.extent.applyTransform(extent, projzh.ll2bmerc),
            tileSize: [256, 256]
          })
        })
      });
    var leibies=['计算机','销售','汽车','贸易/物流','金融','房地产','生产制造','咨询','医疗健康','服务业','管培生/公务员','管理','传媒','其他'];
    var vector_source=[];
    var vectors=[];
    var vectors2=[]
    var color=['#DC143C','#C71585','#4B0082','#000080','#2F4F4F','#228B22','#FFFF00','#800000','#B8860B','#CD853F','#FF0000','#3CB371','#FF1493','#FFA500']
    var num=leibies.length;
    var styleCache = {}
    for (let i=0;i<num;i++){
        var Source = new ol.source.Vector();
        var job_vector=new ol.layer.Vector({
            source:new ol.source.Cluster({
                distance: 30,
                source: Source,
            }),
            style: function (feature, resolution) {
                var size = feature.get('features').length;
                var style = styleCache[size+'_'+i];
                if (!style) {
                    style = [
                        new ol.style.Style({
                            image: new ol.style.Circle({
                                radius: Math.ceil(size+100)/100+5,
                                stroke: new ol.style.Stroke({
                                    color: color[i],
                                    width:2
                                }),
                                fill: new ol.style.Fill({
                                    color: color[i]
                                })
                            }),
                            text: new ol.style.Text({
                                text: size.toString(),
                                fill: new ol.style.Fill({
                                    color: "#262626"
                                })
                            })
                        })
                    ];
                    styleCache[size+'_'+i] = style;
                }
                return style;
            }
        });
        var job_vector2=new ol.layer.Vector({
            source:Source,
            style:new ol.style.Style({
                image: new ol.style.Circle({
                    radius: 1,
                    stroke: new ol.style.Stroke({
                        color: color[i],
                        width:2
                    }),
                    fill: new ol.style.Fill({
                        color: color[i]
                    })
                })
            })

        })
        vector_source.push(Source);
        vectors.push(job_vector);
        vectors2.push(job_vector2);
    }
    var change_num=0;
    var map = new ol.Map({
        target: 'map',
        layers:[baidu],
        view: new ol.View({
            center: ol.proj.fromLonLat([107,39.9299857781]),
            zoom: 5
        })
    });
    var shengfeng=['beijing','shanghai','tianjing','chongqing','guangdong','zhejiang','jiangsu','sichuan','fujian','hainan','shandong',
    'jiangxi','anhui','hebei','henan','hubei','hunan','shanxi','shanxi2','heilongjiang','liaoning','jiling','guangxi','yunnan','gansu',
    'neimenggu','ningxia','xizang','xingjiang','qinghai']     
    var sheng=['北京','上海','天津','重庆','广东','浙江','江苏','四川','福建','海南','山东',
    '江西','安徽','河北','河南','湖北','湖南','陕西','山西','黑龙江','辽宁','吉林','广西','云南','甘肃',
    '内蒙古','宁夏','西藏','新疆','青海']   
    var numsum=new Array(num*shengfeng.length);
    var layers=[]
    map.on('click', showInfo);
    map.on('moveend',change_vector);
    var info = document.getElementById('info');
    function showInfo(event) {
        var features = map.getFeaturesAtPixel(event.pixel);
        if (!features) {
            info.innerText = '';
            info.style.opacity = 0;
            return;
        }
        for(i=0;i<features.length;i++){
        var properties = features[i].getProperties();
        var id=features[i].getId().split(".")[1];
        var url='https://jobs.51job.com/quanguo/'+id+'.html?s=01&t=0';
        var a=document.createElement("a");
        a.href=url;
        a.innerHTML=url;
        a.style.color="orange";
        a.target="_blank";
        info.innerText=properties['job_name']+" \n "+ properties['company']+" \n "+properties['min_salary']+"-"+properties['max_salary']+" \n "+properties['place_name']+" \n "+properties['welfare']+" \n "+properties['time']+" \n ";
        info.style.opacity = 1;
        info.style.left=event.pixel[0]+'px';  
        info.style.top=event.pixel[1]+'px';  
        info.style.position='absolute'; 
        info.appendChild(a);
        }
        }
    function change_vector(){
        if(map.getView().getZoom()>13&& change_num==1){
            for(i=0;i<num;i++){
                map.removeLayer(vectors[i]);
                map.addLayer(vectors2[i]);
            }
            change_num=0
        }
        else if(map.getView().getZoom()<=13&&change_num==0) {
            for(i=0;i<num;i++){
                map.removeLayer(vectors2[i]);
                map.addLayer(vectors[i]);
            }
            change_num=1
        }
    }
    function getprovince(s1){
      var pnum=0;
      for(i=0;i<shengfeng.length;i++){
        if(s1==shengfeng[i]){
          pnum=i;
        }
      }
      return pnum;
    }

    function getjob(s2){
      var jnum=0;
      for(i=0;i<num;i++){
        if(s2==leibies[i]){
          jnum=i;
        }
      }
      return jnum;
    }
    function creat_cqlurl(option1,option2,option3,option4,option5){
        var featureRequest = new ol.format.WFS().writeGetFeature({
            srsName: 'EPSG:3857',
            featureNS: 'http://www.opengeospatial.net/cite',

            featurePrefix: 'cite',
            featureTypes: ['job_point_'+option1],
            outputFormat: 'application/json',
            filter:new ol.format.filter.and(ol.format.filter.equalTo('class2',option2 ),
             ol.format.filter.equalTo('status',true),
             ol.format.filter.or(ol.format.filter.like('class2','*'+option3+'*'),
             ol.format.filter.like('job_name','*'+option3+'*'),
             ol.format.filter.like('place_name','*'+option3+'*'),
             ol.format.filter.like('company','*'+option3+'*')),
             ol.format.filter.and(ol.format.filter.greaterThan("max_salary",parseFloat(option4)),
             ol.format.filter.lessThan("min_salary",parseFloat(option5)))
             )
          }); 
          fetch('http://122.51.15.92:8080/geoserver/wfs', {
            method: 'POST',
            body: new XMLSerializer().serializeToString(featureRequest)
          }).then(function(response) {
            return response.json();
          }).then(function(json) {
            var features = new ol.format.GeoJSON().readFeatures(json);
            var flen=features.length;//这个length只能得到某个省的所有？
            var a=getprovince(option1);//得到行数和列数，省？职业？
            var b=getjob(option2);
            var c=a*num+b;
            numsum[c]=flen;
            for(i=0;i<num;i++){
                if (option2==leibies[i]){
                    vector_source[i].addFeatures(features);
                }
            }
          });
    }
    function search(){
      for(i=0;i<numsum.length;i++){
        numsum[i]=0;
      } 
        var options1= $('#shengfeng').find(".selected");
        var options2=$('#leibie').find(".selected");
        var option3=document.getElementById("g1").value;
        // var option4=document.getElementById("min_s").value;
        // var option5=document.getElementById("max_s").value;
        var option4=0;
        var option5=0;
        var myselect=document.getElementById("sal").value;
        for(i=0;i<num;i++){
            vector_source[i].clear();
        }
        if(options1.length==0){
            options1=$('#shengfeng').find(".fs-option");
        }
        if(options2.length==0){
            options2=$('#leibie').find(".fs-option");
        }
        if(myselect=="salary1"){
          option4='0';
          option5='3000';
        }
        if(myselect=="salary2"){
          option4='3000';
          option5='6000';
        }
        if(myselect=="salary3"){
          option4='6000';
          option5='10000';
        }
        if(myselect=="salary4"){
          option4='10000';
          option5='15000';
        }
        
        if(myselect=="salary5"){
          option4='15000';
          option5='20000';
        }
        //if(op
        if(myselect=="salary6"){
          option4='20000';
          option5='9000000';
        }
        if(myselect=="salary7"){
          option4='0';
          option5='9000000';
        }
        // if(document.getElementById("min_s").value==""){
        //    option4='0';
        // }
        // if(document.getElementById("max_s").value==""){
        //     option5='9000000';
        //  }
        for(i=0;i<options1.length;i++){
            for(n=0;n<options2.length;n++){
                var shengfeng_s=options1[i].getAttribute('data-value');
                var leibie_s=options2[n].getAttribute('data-value');
                if(shengfeng_s!='all'&&leibie_s!='all'){
                    creat_cqlurl(shengfeng_s,leibie_s,option3,option4,option5);
                }
            }
        }  
    }


    // function createpie(){//功能为统计各省加在一起的不同类别占比
    //   var leinumber=[];
    //   for (i=0;i<num;i++){
    //       var features=vector_source[i].getFeatures();
    //       leinumber[i]=features.length;
    //   }
    //   var cal=0;
    //      for(i=0;i<num;i++){
    //          if(leinumber[i]!=0){
    //              cal++;
    //          }
    //      }
    //     var ans1=new Array(cal);
    //     var name1=new Array(cal);

    //     var temp=0;
    //     for(i=0;i<num;i++){
    //         if(leinumber[i]!=0){
    //             ans1[temp]=leinumber[i];
    //             name1[temp]=leibies[i];
    //             temp++;
    //         }
    //     }
    //     var consequense=new Array(cal);
    //     for(i=0;i<cal;i++){
    //         consequense[i]={value:ans1[i], name:name1[i]};
    //     }
    //     var array=[];
    //     for(var i=0;i<cal;i++){
    //       var item={
    //         name:ans1[i],
    //         value:name1[i]
    //       }
    //       array.push(item);
    //     }
    //     //alert(name1[0]);
    //     var option = {
    //     title : {
    //         text: '饼状图',
    //         subtext: '职业类别统计',
    //         x:'center'
    //     },
    //     tooltip : {
    //         trigger: 'item',
    //         formatter: "{b} : {c} ({d}%)"
    //     },
    //     legend: {
    //         orient: 'vertical',
    //         left: 'left',
    //         data: []
    //     },
    //     series : [
    //         {
    //             //name: '职位类别',
    //             type: 'pie',
    //             radius : '55%',
    //             center: ['50%', '60%'],
    //             data:[
    //             ],
    //             itemStyle: {
    //                 emphasis: {
    //                     shadowBlur: 10,
    //                     shadowOffsetX: 0,
    //                     shadowColor: 'rgba(0, 0, 0, 0.5)'
    //                 }
    //             }
    //         }
    //     ]
    // };
    // for(i=0;i<cal;i++){
    //     option.series[0].data[i]=consequense[i];
    //   }
    //   for(i=0;i<cal;i++){
    //     option.legend.data[i]=name1[i];
    //   }
    //     mychart.setOption(option);
    // }

    // function createradar(){//难点在如何与省份匹配
    //   var options1= $('#shengfeng').find(".selected");
    //   var options2= $('#leibie').find(".selected");
    //   var qnmd=options2[0].getAttribute('data-value');
    //   var wc=options1[0].getAttribute('data-value');
    //   var nmsl=getjob(qnmd);
    //   var cnm= getprovince(wc);
    //   var pro=[];
    //   var jo=[];
      
    //   var max1=new Array(options2.length);
    //   var count0=0;
    //   for(j=0;j<num;j++){
    //     //这里错了，根本就没选第一个省份
    //     max1[count0]=0;
    //     if(numsum[cnm*num+j]!=0){
    //       for(i=0;i<shengfeng.length;i++){
    //         if(max1[count0]<numsum[i*num+j]){
    //            max1[count0]=numsum[i*num+j]; 
    //         }
    //       }
    //       count0++;
    //     }
    //   }
    //   for(i=0;i<options2.length;i++){
    //     max1[i]=Math.round(max1[i]*1.15);
    //   }
      
    //   //省份这个字段在哪里提取？
    //   var count1=0;
    //   var count2=0;
    //   for(i=0;i<shengfeng.length;i++){//这里不能是省份的长度
    //     if(numsum[i*num+nmsl]!=0){
    //       pro[count1]=sheng[i];
    //       count1++;
    //     }
    //   }
    //   for(i=0;i<num;i++){
    //     if(numsum[cnm*num+i]!=0){
    //       jo[count2]=leibies[i];
    //       count2++;
    //     }
    //   }
     
    //     option = {
    //       title: {
    //           text: '省份雷达图'
    //       },
    //       tooltip: {},
    //       legend: {
    //           data: []//
    //       },
    //       radar: {
    //           // shape: 'circle',
    //           name: {
    //               textStyle: {
    //                   color: '#fff',
    //                   backgroundColor: '#999',
    //                   borderRadius: 3,
    //                   padding: [3, 5]
    //              }
    //           },
    //           indicator: [
                 
    //           ]
    //       },
    //       series: [{
    //           name: '省份职业对比',
    //           type: 'radar',
    //           // areaStyle: {normal: {}},
    //           data : [
    //               // {
    //               //     value : [4300, 10000, 28000, 35000, 50000, 19000],
    //               //     name : '预算分配（Allocated Budget）'
    //               // },
    //               //  {
    //               //     value : [5000, 14000, 28000, 31000, 42000, 21000],
    //               //     name : '实际开销（Actual Spending）'
    //               // }
    //           ]
    //       }]
    //   };
    //   for(i=0;i<options1.length;i++){
    //     option.legend.data[i]=pro[i];
    //   }
    //   for(i=0;i<options2.length;i++){
    //     option.radar.indicator[i]={name:jo[i],max:max1[i]};
    //   }
    //   var a=0;//表示省份
    //   for(i=0;i<shengfeng.length;i++){
    //     var val=new Array(options2.length);
    //     var b=0;//表示不同职业
    //     for(j=0;j<num;j++){
    //       if(numsum[i*num+j]!=0){//是不是真的是0
    //         val[b]=numsum[i*num+j];
    //         b++;
    //       }
    //     }
    //     option.series[0].data[a]={value:val,name:pro[1]};
    //     //option.legend.data[a]=pro[a];
    //     a++;
    //   }
    //   alert(pro[1]);
    //     mychart.setOption(option);
    // }