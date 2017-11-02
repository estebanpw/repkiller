#pragma once

#include <iostream>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <queue>
#include <thread>
#include <string>

#include "structs.h"
#include "commonFunctions.h"

using namespace std;

class SaverQueue {
private:
  struct SaveRequest {
    SaveRequest(const std::string & path, FGList * fgl) {
      this->path = std::string(path);
      this->fgl = fgl;
    }
    std::string path;
    FGList * fgl;
  };

  mutable size_t count_ = 0;
  bool running_ = false;
  const sequence_manager & seq_mngr;
  std::mutex mutex_;
  std::condition_variable cond_;
  std::queue<SaveRequest> queue_;
  std::unique_ptr<std::thread> thread_ptr_;

  void run();
public:
  SaverQueue(const sequence_manager & seq_mngr) : seq_mngr(seq_mngr) {};
  ~SaverQueue();
  void start();
  void stop();
  void addRequest(const string & path, FGList * fgl);
};
