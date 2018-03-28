#ifndef JSON_PARSER
#define JSON_PARSER

#include <fstream>
#include <string>
#include <exception>

#include "json.hpp"

namespace json_parser{
  
  using json = nlohmann::json;
  
  inline const json& get_child(const json &root, const std::string &name){
    return root[name];
  }
  
  template<typename T>
    inline T get(const json &root, const std::string &name){
    return root[name].get<T>();
  }

  template<typename T>
  inline T get(const json &root, const std::string &name, std::initializer_list<T> acceptable_list){
    T item = root[name].get<T>();
    bool acceptable = false;
    for(auto &dmy : acceptable_list) if(item == dmy) acceptable = true;
    if(!acceptable){
      std::cerr << "Got : " << item << std::endl;
      std::cerr << "Expected one of : ";
      for(auto &dmy : acceptable_list) std::cerr << dmy << ", ";
      std::cerr << std::endl;
      
      std::string buffer = "Invalid value for " + name;
      throw std::domain_error(buffer);
    }
    return item;
  }
  
  template<typename T>
  inline std::vector<T> get_vector(const json &root, const std::string &name,
				   const unsigned long long &size){
    std::vector<T> vec = root[name].get<std::vector<T>>();

    auto dmy_size = vec.size();
    if(dmy_size != size){
      std::string buffer =  "vector " + name + " should be of size " + std::to_string(size)
	+ ". Got " + std::to_string(dmy_size);
      throw std::domain_error(buffer);
    }
    return vec;
  }
  template<typename T>
    inline std::vector<std::vector<T>> get_matrix(const json &root, const std::string &name,
						  const unsigned long long &size1,
						  const unsigned long long &size2){
    std::vector<std::vector<T>> mat = root[name].get<std::vector<std::vector<T>>>();
    if(mat.size() != size1){
      std::string buffer = "matrix " + name + " should have " + std::to_string(size1)
	+ " rows. Got " + std::to_string(mat.size());
      throw std::domain_error(buffer);
    }
    for(auto i = 0; i < mat.size(); i++){
      if(mat[i].size() != size2){
	std::string buffer = "matrix " + name + " should have " + std::to_string(size2)
	  + " columns. Got " + std::to_string(mat.size()) + " for row " + std::to_string(i);
	throw std::domain_error(buffer);
      }
    }
    return mat;
  }
  
  inline void parse_file(const char* fname, json &root){
    std::ifstream input_file(fname);
    root << input_file;
    std::cerr << "# Using JsonCons" << std::endl;    
  }
  inline void parse(const std::string &raw, json &root){
    root = json::parse(raw);
    std::cerr << "# Using JsonCons" << std::endl;
  }
  
  inline void dump(const char* fname, const json &root){
    std::ofstream output_file(fname);
    output_file << root.dump(4) << std::endl;
    output_file.close();
  }
}

#endif
