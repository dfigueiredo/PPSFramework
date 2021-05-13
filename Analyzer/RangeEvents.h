#include <map>
#include <string>

class RangeEvents
{

  private:
    map<std::string, double> mapRange;

    double value0min=0;
    double value1min=0;
    double value2min=0;
    double value3min=0;
    double value4min=0;
    double value5min=0;
    double value6min=0;
    double value7min=0;
    double value8min=0;
    double value9min=0;
    double value10min=0;
    double value11min=0;
    double value12min=0;
    double value13min=0;
    double value14min=0;
    double value15min=0;
    double value16min=0;
    double value17min=0;
    double value18min=0;
    double value0max=0;
    double value1max=0;
    double value2max=0;
    double value3max=0;
    double value4max=0;
    double value5max=0;
    double value6max=0;
    double value7max=0;
    double value8max=0;
    double value9max=0;
    double value10max=0;
    double value11max=0;
    double value12max=0;
    double value13max=0;
    double value14max=0;
    double value15max=0;
    double value16max=0;
    double value17max=0;
    double value18max=0;
    double value19max=0;
    double entries = 0;
    std::string key;

    double getMap(std::string);

  public:
    RangeEvents(std::string key_, double entries_){
      key = key_;
      entries = entries_;
      value0min = getMap(key+"_0");
      value1min = value0min + getMap(key+"_1");
      value2min = value1min + getMap(key+"_2");
      value3min = value2min + getMap(key+"_3");
      value4min = value3min + getMap(key+"_4");
      value5min = value4min + getMap(key+"_5");
      value6min = value5min + getMap(key+"_6");
      value7min = value6min + getMap(key+"_7");
      value8min = value7min + getMap(key+"_8");
      value9min = value8min + getMap(key+"_9");
      value10min = value9min + getMap(key+"_10");
      value11min = value10min + getMap(key+"_11");
      value12min = value11min + getMap(key+"_12");
      value13min = value12min + getMap(key+"_13");
      value14min = value13min + getMap(key+"_14");
      value15min = value14min + getMap(key+"_15");
      value16min = value15min + getMap(key+"_16");
      value17min = value16min + getMap(key+"_17");
      value18min = value17min + getMap(key+"_18");
      value0max = floor(entries*getMap(key+"_0"));
      value1max = floor(entries*value0min)+floor(entries*getMap(key+"_1"));
      value2max = floor(entries*value1min)+floor(entries*getMap(key+"_2"));
      value3max = floor(entries*value2min)+floor(entries*getMap(key+"_3"));
      value4max = floor(entries*value3min)+floor(entries*getMap(key+"_4"));
      value5max = floor(entries*value4min)+floor(entries*getMap(key+"_5"));
      value6max = floor(entries*value5min)+floor(entries*getMap(key+"_6"));
      value7max = floor(entries*value6min)+floor(entries*getMap(key+"_7"));
      value8max = floor(entries*value7min)+floor(entries*getMap(key+"_8"));
      value9max = floor(entries*value8min)+floor(entries*getMap(key+"_9"));
      value10max = floor(entries*value9min)+floor(entries*getMap(key+"_10"));
      value11max = floor(entries*value10min)+floor(entries*getMap(key+"_11"));
      value12max = floor(entries*value11min)+floor(entries*getMap(key+"_12"));
      value13max = floor(entries*value12min)+floor(entries*getMap(key+"_13"));
      value14max = floor(entries*value13min)+floor(entries*getMap(key+"_14"));
      value15max = floor(entries*value14min)+floor(entries*getMap(key+"_15"));
      value16max = floor(entries*value15min)+floor(entries*getMap(key+"_16"));
      value17max = floor(entries*value16min)+floor(entries*getMap(key+"_17"));
      value18max = floor(entries*value17min)+floor(entries*getMap(key+"_18"));
      value19max = floor(entries*value18min)+floor(entries*getMap(key+"_19"));
    };
    virtual ~RangeEvents() {}
    double getMinimum(int, std::string);
    double getMaximum(int, std::string);

};

double RangeEvents::getMap(std::string key) 
{

  double value = 0;

  mapRange["muon_0"]  = 0.009941;   //0.0135092, Fractions for 1 Multi-RP and 1 pimapRangeels
  mapRange["muon_1"]  = 0.011229;   //0.0150113;
  mapRange["muon_2"]  = 0.016891;   //0.0222738;
  mapRange["muon_3"]  = 0.018763;   //0.0244954;
  mapRange["muon_4"]  = 0.065784;   //0.0847967;
  mapRange["muon_5"]  = 0.097749;   //0.119802;
  mapRange["muon_6"]  = 0.064967;   //0.0747702;
  mapRange["muon_7"]  = 0.093342;   //0.103186;
  mapRange["muon_8"]  = 0.019649;   //0.021513;
  mapRange["muon_9"]  = 0.041326;   //0.0417492;
  mapRange["muon_10"] = 0.041221;   //0.0389873;
  mapRange["muon_11"] = 0.042031;   //0.0375865;
  mapRange["muon_12"] = 0.023791;   //0.0297421;
  mapRange["muon_13"] = 0.019781;   //0.023547;
  mapRange["muon_14"] = 0.000317;   //0.000372487;
  mapRange["muon_15"] = 0.008085;   //0.00945544;
  mapRange["muon_16"] = 0.185997;   //0.174744;
  mapRange["muon_17"] = 0.098481;   //0.0764684;
  mapRange["muon_18"] = 0.013582;   //0.0125681;
  mapRange["muon_19"] = 0.104090;   //0.0754216;

  mapRange["electron_0"]  = 0.015610;   //0.0204228;     
  mapRange["electron_1"]  = 0.022636;   //0.0292415;
  mapRange["electron_2"]  = 0.031372;   //0.0399177;
  mapRange["electron_3"]  = 0.039083;   //0.0491031;
  mapRange["electron_4"]  = 0.072814;   //0.090098;
  mapRange["electron_5"]  = 0.119773;   //0.140086;
  mapRange["electron_6"]  = 0.090694;   //0.0986203;
  mapRange["electron_7"]  = 0.144755;   //0.151347;
  mapRange["electron_8"]  = 0.013952;   //0.014786;
  mapRange["electron_9"]  = 0.031176;   //0.0300992;
  mapRange["electron_10"] = 0.032689;   //0.0298066;
  mapRange["electron_11"] = 0.038521;   //0.0299309;
  mapRange["electron_12"] = 0.017496;   //0.0211467;
  mapRange["electron_13"] = 0.015451;   //0.0176158;
  mapRange["electron_14"] = 0.000225;   //0.000250771;
  mapRange["electron_15"] = 0.005792;   //0.00657242;
  mapRange["electron_16"] = 0.131171;   //0.117324;
  mapRange["electron_17"] = 0.072942;   //0.0537374;
  mapRange["electron_18"] = 0.009866;   //0.00885834;
  mapRange["electron_19"] = 0.074057;   //0.0510357;

  value = mapRange.find(key)->second;
  return(value); 

}

double RangeEvents::getMinimum(int xangle, std::string era){

  double rangemin = 0; 

  if(xangle==120 && era=="B"){
    rangemin=0;
  }
  if(xangle==130 && era=="B"){
    rangemin=floor(entries*value0min);
  }
  if(xangle==140 && era=="B"){
    rangemin=floor(entries*value1min);
  }
  if(xangle==150 && era=="B"){
    rangemin=floor(entries*value2min);
  }
  if(xangle==120 && era=="C"){
    rangemin=floor(entries*value3min);
  }
  if(xangle==130 && era=="C"){
    rangemin=floor(entries*value4min);
  }
  if(xangle==140 && era=="C"){
    rangemin=floor(entries*value5min);
  }
  if(xangle==150 && era=="C"){
    rangemin=floor(entries*value6min);
  }
  if(xangle==120 && era=="D"){
    rangemin=floor(entries*value7min);
  }
  if(xangle==130 && era=="D"){
    rangemin=floor(entries*value8min);
  }
  if(xangle==140 && era=="D"){
    rangemin=floor(entries*value9min);
  }
  if(xangle==150 && era=="D"){
    rangemin=floor(entries*value10min);
  }
  if(xangle==120 && era=="E"){
    rangemin=floor(entries*value11min);
  }
  if(xangle==130 && era=="E"){
    rangemin=floor(entries*value12min);
  }
  if(xangle==140 && era=="E"){
    rangemin=floor(entries*value13min);
  }
  if(xangle==150 && era=="E"){
    rangemin=floor(entries*value14min);
  }
  if(xangle==120 && era=="F"){
    rangemin=floor(entries*value15min);
  }
  if(xangle==130 && era=="F"){
    rangemin=floor(entries*value16min);
  }
  if(xangle==140 && era=="F"){
    rangemin=floor(entries*value17min);
  }
  if(xangle==150 && era=="F"){
    rangemin=floor(entries*value18min);
  }

  return(rangemin);

}

double RangeEvents::getMaximum(int xangle, std::string era){

  double rangemax = 0;

  if(xangle==120 && era=="B"){
    rangemax=value0max;
  }
  if(xangle==130 && era=="B"){
    rangemax=value1max;
  }
  if(xangle==140 && era=="B"){
    rangemax=value2max;
  }
  if(xangle==150 && era=="B"){
    rangemax=value3max;
  }
  if(xangle==120 && era=="C"){
    rangemax=value4max;
  }
  if(xangle==130 && era=="C"){
    rangemax=value5max;
  }
  if(xangle==140 && era=="C"){
    rangemax=value6max;
  }
  if(xangle==150 && era=="C"){
    rangemax=value7max;
  }
  if(xangle==120 && era=="D"){
    rangemax=value8max;
  }
  if(xangle==130 && era=="D"){
    rangemax=value9max;
  }
  if(xangle==140 && era=="D"){
    rangemax=value10max;
  }
  if(xangle==150 && era=="D"){
    rangemax=value11max;
  }
  if(xangle==120 && era=="E"){
    rangemax=value12max;
  }
  if(xangle==130 && era=="E"){
    rangemax=value13max;
  }
  if(xangle==140 && era=="E"){
    rangemax=value14max;
  }
  if(xangle==150 && era=="E"){
    rangemax=value15max;
  }
  if(xangle==120 && era=="F"){
    rangemax=value16max;
  }
  if(xangle==130 && era=="F"){
    rangemax=value17max;
  }
  if(xangle==140 && era=="F"){
    rangemax=value18max;
  }
  if(xangle==150 && era=="F"){
    rangemax=value19max;
  }

  return(rangemax);

}
