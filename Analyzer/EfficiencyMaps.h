#include <map>
#include <string>

class EfficiencyMaps
{
  private:
    // Internally we use a map.
    // But typedef the class so it is easy to refer too.
    // Also if you change the type you only need to do it in one place.
    typedef std::map<std::string, std::string>  effMap;
    effMap   x; // The data store.

    // The only copy of the map
    // I dont see any way of modifying so declare it const (unless you want to modify it)
    static const Mapping myMap;

    // Make the constructor private.
    // This class is going to hold the only copy.
    EfficiencyMaps()
    {
      x["muon_0"]  = 0.009941;   //0.0135092, Fractions for 1 Multi-RP and 1 pixels
      x["muon_1"]  = 0.011229;   //0.0150113;
      x["muon_2"]  = 0.016891;   //0.0222738;
      x["muon_3"]  = 0.018763;   //0.0244954;
      x["muon_4"]  = 0.065784;   //0.0847967;
      x["muon_5"]  = 0.097749;   //0.119802;
      x["muon_6"]  = 0.064967;   //0.0747702;
      x["muon_7"]  = 0.093342;   //0.103186;
      x["muon_8"]  = 0.019649;   //0.021513;
      x["muon_9"]  = 0.041326;   //0.0417492;
      x["muon_10"] = 0.041221;   //0.0389873;
      x["muon_11"] = 0.042031;   //0.0375865;
      x["muon_12"] = 0.023791;   //0.0297421;
      x["muon_13"] = 0.019781;   //0.023547;
      x["muon_14"] = 0.000317;   //0.000372487;
      x["muon_15"] = 0.008085;   //0.00945544;
      x["muon_16"] = 0.185997;   //0.174744;
      x["muon_17"] = 0.098481;   //0.0764684;
      x["muon_18"] = 0.013582;   //0.0125681;
      x["muon_19"] = 0.104090;   //0.0754216;

      x["electron_0"]  = 0.015610;   //0.0204228;     
      x["electron_1"]  = 0.022636;   //0.0292415;
      x["electron_2"]  = 0.031372;   //0.0399177;
      x["electron_3"]  = 0.039083;   //0.0491031;
      x["electron_4"]  = 0.072814;   //0.090098;
      x["electron_5"]  = 0.119773;   //0.140086;
      x["electron_6"]  = 0.090694;   //0.0986203;
      x["electron_7"]  = 0.144755;   //0.151347;
      x["electron_8"]  = 0.013952;   //0.014786;
      x["electron_9"]  = 0.031176;   //0.0300992;
      x["electron_10"] = 0.032689;   //0.0298066;
      x["electron_11"] = 0.038521;   //0.0299309;
      x["electron_12"] = 0.017496;   //0.0211467;
      x["electron_13"] = 0.015451;   //0.0176158;
      x["electron_14"] = 0.000225;   //0.000250771;
      x["electron_15"] = 0.005792;   //0.00657242;
      x["electron_16"] = 0.131171;   //0.117324;
      x["electron_17"] = 0.072942;   //0.0537374;
      x["electron_18"] = 0.009866;   //0.00885834;
      x["electron_19"] = 0.074057;   //0.0510357;
 
    }


  public:
    // Public interface.
    //    Returns a const reference to the value.
    //    The interface use static methods (means we dont need an instance)
    //    Internally we refer to the only instance.
    static std::string const& getValue(std::string const& value)
    {
      // Use find rather than operator[].
      // This way you dont go inserting garbage into your data store.
      // Also it allows the data store to be const (as operator may modify the data store
      // if the value is not found).

      effMap::const_iterator   find    = myMap.x.find(value);
      if (find != myMap.x.end())
      {
	// If we find it return the value.
	return find->second;
      }

      // What happens when we don;t find anything.
      // Your original code created a garbage entry and returned that.
      // Could throw an exception or return a temporary reference.
      // Maybe ->  throw int(1);
      return "";
    }

};
