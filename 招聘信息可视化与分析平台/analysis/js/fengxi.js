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
var shengfeng=['beijing','shanghai','tianjing','chongqing','guangdong','zhejiang','jiangsu','sichuan','fujian','hainan','shandong',
    'jiangxi','anhui','hebei','henan','hubei','hunan','shanxi','shanxi2','heilongjiang','liaoning','jiling','guangxi','yunnan','gansu',
    'neimenggu','ningxia','xizang','xingjiang','qinghai'];
var shengfeng2=['北京','上海','天津','重庆','广东','浙江','江苏','四川','福建','海南','山东',
    '江西','安徽','河北','河南','湖北','湖南','陕西','山西','黑龙江','辽宁','吉林','广西','云南','甘肃',
    '内蒙古','宁夏','西藏','新疆','青海'];

var all_source=new ol.source.Vector();
var heatmap=new ol.layer.Heatmap({
    source:new ol.source.Cluster({
        distance:5,
        source: all_source,
    }),
    blur: 30,
    radius: 10
});
var map = new ol.Map({
    target: 'map',
    layers:[baidu,heatmap],
    view: new ol.View({
        center: ol.proj.fromLonLat([113.734427,23.210562]),
        zoom: 5
    })
});
var dom1 = document.getElementById("c1");
var dom2 = document.getElementById("c2");
var dom3 = document.getElementById("c3");
var myChart1 = echarts.init(dom1);
var myChart2 = echarts.init(dom2);
var myChart3 = echarts.init(dom3);
var data=[];
var jingdu=0;
var features={}

function creat_cqlurl(times){
  function dataFormatter(obj) {
    var temp;
    for (let time = 0; time <times.length; time++) {
        var sum = 0;
        temp = obj[times[time]];
        for (var i = 0, l = temp.length; i < l; i++) {
            sum += temp[i];
            obj[times[time]][i] = {
                name : shengfeng2[i],
                value : temp[i]
            }
        }
        obj[times[time] + 'sum'] = sum;
    }
    return obj;
  }
  data=[];
  all_source.clear()

  for(let i=0;i<shengfeng.length;i++){
    for(let idx=0;idx<times.length;idx++){
      var featureRequest = new ol.format.WFS().writeGetFeature({
        srsName: 'EPSG:3857',
        featureNS: 'http://www.opengeospatial.net/cite',
        featurePrefix: 'cite',
        featureTypes: [shengfeng[i]+'_view'],
        outputFormat: 'application/json',
        filter:ol.format.filter.equalTo('time',times[idx]),
      }); 
        fetch('http://122.51.15.92:8080/geoserver/wfs', {
        method: 'POST',
        body: new XMLSerializer().serializeToString(featureRequest)
      }).then(function(response) {
        return response.json();
      }).then(function(json) {
        var features = new ol.format.GeoJSON().readFeatures(json);
        if(features.length>0){
        data.push({shengfeng:shengfeng[i],time:times[idx],feature:features});
        }
        jingdu--;
        $('.progress-bar').css('width',parseInt((1-jingdu/(shengfeng.length*times.length))*100)+'%');
        var observer = new MutationObserver(function(mutations) {
          mutations.forEach(function(mutation) {
            if (mutation.type == 'attributes' && mutation.attributeName == 'style') {
              var el = mutation.target;
                var width = el.style.width; // Can't use jQuery here, as we need the value back in percent
            var $parentEl =$(el).parent('.progress');
            $parentEl.attr('data-width',width); // Why doesn't this work?? $parentEl.data('width',width)
            $parentEl.find('.progress-text').text(width);
            }
          });
        });
        var config = {
          attributes: true,
          attributeFilter: ['style'],
          childList: false,
          characterData: false
        };
      
        $('.progress-bar').each(function(e) {
          observer.observe(this, config);
        })
      
        if(jingdu==0){
          var dvObj = document.getElementById('jz');
          dvObj.innerHTML = ''
          var c1_data=[];
          var c2_data={};
          var c3_data={};
          var time3=times[0].split('-');
          var time4=times[times.length-1].split('-');
          var timeplay = new TimePlay({
            speed: 2000,//播放速度
            startDate: parseInt(time3[0]+time3[1]+time3[2])-1,//开始日
            endDate: parseInt(time4[0]+time4[1]+time4[2]),//结束日期
            timeUnitControl: false,//是否显示时/天切换控件
            onClickChangeEnd: function(){//点击后回调
              //获取点击的时间
              var day  = timeplay.curr_day.day,//日
              mon  = timeplay.curr_day.mon,//月
              year = timeplay.curr_day.year;//年
              console.log(year+"年" + mon + "月" + day + "日" );
            },
            onAnimateEnd: function(){//时间轴动画每次结束回调
              var day  = timeplay.curr_day.day,//日
              mon  = timeplay.curr_day.mon,//月
              year = timeplay.curr_day.year;//年
              console.log(year+"年" + mon + "月" + day + "日" );
            }
          });
          for(let i1=0;i1<shengfeng2.length;i1++){
            c3_data[shengfeng2[i1]]=[];
            for(let i2=0;i2<times.length;i2++){
              c3_data[shengfeng2[i1]].push(0);
            }
          }
          for(let i1=0;i1<leibies.length;i1++){
            c1_data.push([])
            c2_data[leibies[i1]]={};
            for(let i2=0;i2<shengfeng.length;i2++){
              c1_data[i1].push(0)
            }
            for(let i3=0;i3<times.length;i3++){
              c2_data[leibies[i1]][times[i3]]=[];
              for(let i4=0;i4<shengfeng.length;i4++){
                c2_data[leibies[i1]][times[i3]].push(0)
              }
            }
            }
          
          for(z=0;z<data.length;z++){
            if(data[z]['time']==times[0])
            { 
              all_source.addFeatures(data[z]['feature']);
          }  
            for(let i1=0;i1<data[z]['feature'].length;i1++){
            for(let i2=0;i2<shengfeng.length;i2++){
              for(let i3=0;i3<leibies.length;i3++){
                  if(data[z]['feature'][i1].getProperties()['class2']==leibies[i3]&&data[z]['shengfeng']==shengfeng[i2]){
                    c1_data[i3][i2]+=1;
                    c2_data[leibies[i3]][data[z]['time']][i2]+=1;
                    for(let i4=0;i4<times.length;i4++){
                      if(data[z]['time']==times[i4]){
                        c3_data[shengfeng2[i2]][i4]+=1;
                        break;
                      }
                    }

                    continue;
                  }
              }
              continue;
            }
            }
        } 
        var start_times=times[0].split('-');
        var start_time='';
        var end_times=times[times.length-1].split('-');
        var end_time='';
        for(let i1=0;i1<3;i1++){
          start_time+=start_times[i1];
          end_time+=end_times[i1];
        }
        $("#pause").click(function(){
          timeplay.delayAnimation();//延迟动画
        })
        
        $("#play").click(function(){
          timeplay.continueAnimation();//继续动画
        })
        for(let i5=0;i5<leibies.length;i5++){
          c2_data[leibies[i5]]=dataFormatter(c2_data[leibies[i5]]);
        }
        var chart_d=[];
        for(let i1=0;i1<shengfeng.length;i1++){
        var childers=[];
        for(let i2=0;i2<leibies.length;i2++){
          childers.push({name:leibies[i2]+c1_data[i2][i1],value:c1_data[i2][i1],itemStyle:{color:'#3e0317'}})
        }
        chart_d.push({name:shengfeng2[i1],itemStyle:{color:'#f99e1c'},children:childers})
          }
          var option={
          series: {
          type: 'sunburst',
          highlightPolicy: 'ancestor',
          data: chart_d,
          radius: [0, '90%'],
          label: {
              rotate: 'radial',
          }
      }
          }
        var option2=[]
        for(let i6=0;i6<times.length;i6++){
              var series=[];
              var child=[];
              for(let i7=0;i7<leibies.length;i7++){
                series.push({data:c2_data[leibies[i7]][times[i6]]})
                child.push({name:leibies[i7],value:c2_data[leibies[i7]][times[i6]+'sum']})
              }
              series.push({data:child});
              option2.push({title:{text:times[i6]+'招聘信息数据'},series:series});
        }
        var series2=[];
        for(let i9=0;i9<leibies.length;i9++){
          series2.push({name:leibies[i9],type:'bar'})
        }
        series2.push({
                name: 'GDP占比',
                type: 'pie',
                center: ['75%', '35%'],
                radius: '28%',
                z: 100
            })
        var options2={
          baseOption: {
            timeline: {
                axisType: 'category',
                // realtime: false,
                // loop: false,
                autoPlay: true,
                // currentIndex: 2,
                playInterval: 1000,
                // controlStyle: {
                //     position: 'left'
                // },
                data:times,
                label: {
                }
            },
            title: {
                subtext: '数据来自国家统计局'
            },
            tooltip: {
            },
            legend: {
                x: 'right',
                data: leibies,
            },
            calculable : true,
            grid: {
                top: 80,
                bottom: 100,
                tooltip: {
                    trigger: 'axis',
                    axisPointer: {
                        type: 'shadow',
                        label: {
                            show: true,
                            formatter: function (params) {
                                return params.value.replace('\n', '');
                            }
                        }
                    }
                }
            },
            xAxis: [
                {
                    'type':'category',
                    'axisLabel':{'interval':0},
                    'data':['北京','\n上海','天津','\n重庆','广东','\n浙江','江苏','\n四川','福建','\n海南','山东',
                    '\n江西','安徽','\n河北','河南','\n湖北','湖南','\n陕西','山西','\n黑龙江','辽宁','\n吉林','广西','\n云南','甘肃',
                    '\n内蒙古','宁夏','\n西藏','新疆','\n青海'],
                    splitLine: {show: false}
                }
            ],
            yAxis: [
                {
                    type: 'value',
                    name: 'GDP（亿元）'
                }
            ],
            series:series2
        },
        options:option2
        }
        var series3=[];
        for(let i9=0;i9<shengfeng2.length;i9++){
          series3.push({name:shengfeng2[i9],type:'line',smooth:true,data:c3_data[shengfeng2[i9]]})
        }
        var option3= {
          tooltip: {
              trigger: 'axis',
            
          },
          legend: {
              data:shengfeng2
          },
          grid: {
              top: 70,
              bottom: 50,
              containLabel: true
          },
          toolbox: {
            feature: {
                saveAsImage: {}
            }
        },
          xAxis: 
              {
                  type: 'category',
                  boundaryGap: false,
                  data:times
              },
          yAxis: 
              {
                  type: 'value'
              },
          series:series3
      };
      
          myChart1.setOption(option, true);
          myChart2.setOption(options2, true);
          myChart3.setOption(option3, true);

        }
        })
    }
  }
}
function getBetweenDateStr(start,end){
  var result = [];
  var beginDay = start.split("-");
  var endDay = end.split("-");
  var diffDay = new Date();
  var dateList = new Array;
  var i = 0;
  diffDay.setDate(beginDay[2]);
  diffDay.setMonth(beginDay[1]-1);
  diffDay.setFullYear(beginDay[0]);
  result.push(start);
  while(i == 0&&start!=end){
      var countDay = diffDay.getTime() + 24 * 60 * 60 * 1000;
      diffDay.setTime(countDay);
      dateList[2] = diffDay.getDate();
      dateList[1] = diffDay.getMonth() + 1;
      dateList[0] = diffDay.getFullYear();
      if(String(dateList[1]).length == 1){dateList[1] = "0"+dateList[1]};
      if(String(dateList[2]).length == 1){dateList[2] = "0"+dateList[2]};
      result.push(dateList[0]+"-"+dateList[1]+"-"+dateList[2]);
      if(dateList[0] == endDay[0] && dateList[1] == endDay[1] && dateList[2] == endDay[2]){ i = 1;
      }
  };
  console.log(result);
  return result;
};
function search(){
  var dvObj = document.getElementById('jz');
  dvObj.innerHTML = '<div class="progress aqua" data-width="0%"><div class="progress-text">0%</div><div class="progress-bar"><div class="progress-text">0%</div></div></div>';
  var stime=document.getElementById("date").value;
  var etime=document.getElementById("date2").value;
  var stimes=stime.split('-');
  var etimes=etime.split('-');
  if(parseInt(stimes[1])<10){
    stimes[1]='0'+stimes[1];
  }
  if(parseInt(stimes[2])<10){
    stimes[2]='0'+stimes[2];
  }
  if(parseInt(etimes[1])<10){
    etimes[1]='0'+etimes[1];
  }
  if(parseInt(etimes[2])<10){
    etimes[2]='0'+etimes[2];
  }
  var stime2=stimes[0]+'-'+stimes[1]+'-'+stimes[2];
  var etime2=etimes[0]+'-'+etimes[1]+'-'+etimes[2];
  var times=getBetweenDateStr(stime2,etime2);
  jingdu+=shengfeng.length*times.length;
  creat_cqlurl(times);   
}
