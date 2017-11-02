#include "SaverQueue.h"


void SaverQueue::run() {
  while (true) {
    std::unique_lock<std::mutex> lck(mutex_);
    if (queue_.empty()) {
      if (running_) cond_.wait(lck);
      else return;
    }
    SaveRequest sr = queue_.front();
    queue_.pop();
    lck.unlock();
    try {
      save_all_frag_pairs(sr.path, seq_mngr, *sr.fgl);
    } catch (const runtime_error & e) {
      auto default_path = "represults-" + std::to_string(++count_) + ".csv";
      std::cerr << "Couldn't access " << sr.path << ", saving into " << default_path << "\n" << std::flush;
      save_all_frag_pairs(default_path, seq_mngr, *sr.fgl);
    }
    delete sr.fgl;
  }
}


SaverQueue::~SaverQueue() {
  if (running_) stop();
}

void SaverQueue::start() {
  if (!running_) {
    running_ = true;
    thread_ptr_ = std::unique_ptr<std::thread>(new thread(&SaverQueue::run, this));
  }
}


void SaverQueue::stop() {
  if (running_) {
    running_ = false;
    cond_.notify_all();
    thread_ptr_->join();
  }
}


void SaverQueue::addRequest(const string & path, FGList * fgl) {
  std::unique_lock<std::mutex> lck(mutex_);
  queue_.emplace(path, fgl);
  cond_.notify_all();
}
