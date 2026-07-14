#include "ThreadPool.hpp"
#include <iostream>

ThreadPool::ThreadPool(size_t num_threads_param, ThreadInitFunc init_func) : num_threads(num_threads_param == 0 ? std::thread::hardware_concurrency() : num_threads_param), thread_init(std::move(init_func))
{
    for (size_t i = 0; i < num_threads; ++i) {
        workers.emplace_back([this, i] {
            if (thread_init)
                thread_init(i);
            worker_thread();
        });
    }
}

ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    cv.notify_all();
    
    for (std::thread& worker : workers) {
        if (worker.joinable())
            worker.join();
    }
}

void ThreadPool::worker_thread() {
    while (true) {
        std::unique_lock<std::mutex> lock(queue_mutex);
        cv.wait(lock, [this] { return stop || !tasks.empty(); });
        
        if (stop && tasks.empty())
            return;
        
        if (tasks.empty())
            continue;
        
        auto task = std::move(tasks.front());
        tasks.pop();
        lock.unlock();
        
        try {
            task();
        } catch (const std::exception& e) {
            std::cerr << "Task exception: " << e.what() << std::endl;
        }
        
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            pending_tasks--;
        }
        cv_done.notify_all();
    }
}

void ThreadPool::wait() {
    std::unique_lock<std::mutex> lock(queue_mutex);
    cv_done.wait(lock, [this] { return pending_tasks == 0; });
}
