#pragma once

#include <string>
#include <memory>
#include <initializer_list>
#include <map>
#include <assert.h>

#include <RooRealVar.h>
#include <RooDataSet.h>

class BaseFunction {

  public:
    BaseFunction(const std::string& name):
      mName(name), mHasDataset(false), mHasParameters(false), mIsExtended(false)
  {
  }

    const RooAbsPdf& getPdf() {
      return *mPdf.get();
    }

    virtual bool createPdf(RooRealVar& observable) = 0;

    void setParameters(const std::map<std::string, std::shared_ptr<RooRealVar>>& parameters) {
      mParameters = parameters;
      mHasParameters = true;
    }

    void setDataset(RooDataSet* dataset) {
      mDataset.reset(dataset);
      mHasDataset = true;
    }

    std::map<std::string, std::shared_ptr<RooRealVar>>& getParameters() {
      return mParameters;
    }

    bool isExtended() const {
      return mIsExtended;
    }

  protected:
    std::string mName;
    std::shared_ptr<RooAbsPdf> mPdf;

    std::map<std::string, std::shared_ptr<RooRealVar>> mParameters;
    std::shared_ptr<RooDataSet> mDataset;

    bool mHasDataset;
    bool mHasParameters;
    bool mIsExtended;
};

float jsonValueToFloat(const RooRealVar& observable, int mass, const Json::Value& value) {
  if (value.isConvertibleTo(Json::ValueType::realValue))
    return value.asFloat();

  if (! value.isString()) {
    assert(false);
    return 0.;
  }

  const std::string s = value.asString();

  if (s == "%mass%")
    return mass;

  if (s == "%low_bound%")
    return observable.getBinning().lowBound();

  if (s == "%high_bound%")
    return observable.getBinning().highBound();

  assert(false);
  return 0;
}

std::map<std::string, std::shared_ptr<RooRealVar>> jsonToRealVars(const RooRealVar& observable, int mass, const std::string& prefix, const Json::Value& values) {
  std::map<std::string, std::shared_ptr<RooRealVar>> result;
  std::shared_ptr<RooRealVar> var;

  // Iterator over values for each member
  const Json::Value::Members& names = values.getMemberNames();
  for (const std::string& name: names) {
    if (! values[name].isArray())
      continue;

    const Json::Value& value = values[name];

    std::string title = prefix + "_" + name;

    if (value.size() == 1) {
      var.reset(new RooRealVar(title.c_str(), title.c_str(), jsonValueToFloat(observable, mass, value[0])));
    } else if (value.size() == 3) {
      var.reset(new RooRealVar(title.c_str(), title.c_str(), jsonValueToFloat(observable, mass, value[0]), jsonValueToFloat(observable, mass, value[1]), jsonValueToFloat(observable, mass, value[2])));
    } else {
      assert(false);
    }

    result[name] = var;
  }

  return result;
}
