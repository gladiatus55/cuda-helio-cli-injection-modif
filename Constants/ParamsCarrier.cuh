#ifndef COSMIC_SINGLETON_H
#define COSMIC_SINGLETON_H

#include <string>
#include <sstream>
#include <algorithm>
#include <map>
#include <stdio.h>

class ParamsCarrier{
	struct any {
		enum any_type :char { string_t = 0, int_t = 1, float_t = 2, double_t = 3 };
		any() {
		}
		any(const any& a) {
			this->type = a.type;
			switch (this->type) {
			case any_type::string_t: new(&(this->str)) std::string(a.str); break;
			case any_type::int_t: this->i = a.i; break;
			case any_type::float_t: this->f = a.f; break;
			case any_type::double_t: this->d = a.d; break;
			}
		}
		~any() {
			switch (this->type) {
			case any_type::string_t: { if (str.size()) { str.std::string::~string(); } } break;
			default:;
			}
		}

		any_type type;
		union {
			std::string str;
			int i;
			float f;
			double d;
		};
	};
	using any_t = any::any_type;

private:
	static ParamsCarrier *INSTANCE;
	ParamsCarrier();
	std::map<std::string, any> m;

public:
	static ParamsCarrier *instance(){
		if (!INSTANCE)
			INSTANCE = new ParamsCarrier();
		return INSTANCE;
	}

	void putString(std::string key, std::string value) {
		any a;
		a.type = any_t::string_t;
		new(&(a.str)) std::string(value);
		m.insert({ key ,a });
	}

	void putFloat(std::string key, float value) {
		any a;
		a.type = any_t::float_t;
		a.f = value; 
		m.insert({ key,a });
	}

	void putInt(std::string key, int value) {
		any a;
		a.type = any_t::int_t;
		a.i = value;
		m.insert({ key,a });
	}

	void putDouble(std::string key, double value) {
		any a;
		a.type = any_t::double_t;
		a.d = value;
		m.insert({ key,a });
	}

	int getInt(std::string key, int defaultValue) {
		auto search = m.find(key);
		if (search != m.end()) {
			if (search->second.type == any_t::int_t) {
				return search->second.i;
			}
			return defaultValue; 
		}
		return defaultValue; 
	}

	float getFloat(std::string key, float defaultValue) {
		auto search = m.find(key);
		if (search != m.end()) {
			if (search->second.type == any_t::float_t) {
				return search->second.f;
			}
			return defaultValue;
		}
		return defaultValue;
	}

	double getDouble(std::string key, double defaultValue) {
		auto search = m.find(key);
		if (search != m.end()) {
			if (search->second.type == any_t::double_t) {
				return search->second.d;
			}
			return defaultValue;
		}
		return defaultValue;
	}

	std::string getString(std::string key, std::string defaultValue) {
		auto search = m.find(key);
		if (search != m.end()) {
			if (search->second.type == any_t::string_t) {
				return search->second.str;
			}
			return defaultValue;
		}
		return defaultValue;
	}

};

#endif