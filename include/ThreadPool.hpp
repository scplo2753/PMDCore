#pragma once
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <functional>
#include <memory>

class ThreadPool {
public:
    using ThreadInitFunc = std::function<void(size_t)>;

    explicit ThreadPool(size_t num_threads = 0, ThreadInitFunc init_func = nullptr);
    ~ThreadPool();
    
    // 禁止复制
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    
    template <typename F, typename... Args>
    void enqueue(F&& func, Args&&... args) {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            if (stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");
            
            tasks.emplace(std::bind(std::forward<F>(func), std::forward<Args>(args)...));
            pending_tasks++;
        }
        cv.notify_one();
    }
    
    // 等待所有任务完成
    void wait();
    
    size_t get_thread_count() const { return num_threads; }
    
private:
    void worker_thread();
    
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    
    std::mutex queue_mutex;
    std::condition_variable cv;
    std::condition_variable cv_done;
    
    bool stop = false;
    size_t num_threads = 0;
    size_t pending_tasks = 0;
    ThreadInitFunc thread_init = nullptr;
};
