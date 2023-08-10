#pragma once

#include <mutex>
#include <exception>

namespace pn {

    // See https://stackoverflow.com/questions/11828539/elegant-exception-handling-in-openmp
    class ThreadException {
        std::exception_ptr ptr;
        std::mutex lock;
    public:
        ThreadException() : ptr(nullptr) {}

        ~ThreadException() {
            this->rethrow();
        }

        void rethrow() {
            if (this->ptr) {
                std::exception_ptr e = this->ptr;
                this->ptr = nullptr;
                std::rethrow_exception(e);
            }
        }

        bool hasException() {
            return (bool) this->ptr;
        }

        void captureException() {
            std::unique_lock<std::mutex> guard(this->lock);
            this->ptr = std::current_exception();
        }

        template<typename Function, typename... Parameters>
        void run(Function f, Parameters... params) {
            try {
                f(params...);
            } catch (...) {
                captureException();
            }
        }
    };

#define OMP_DISPATCHER      pn::ThreadException threadExceptionDispatcher;
#define OMP_TRY             try { if (threadExceptionDispatcher.hasException()) continue;
#define OMP_CATCH           } catch (...) { threadExceptionDispatcher.captureException(); }
#define OMP_RETHROW         threadExceptionDispatcher.rethrow();


}
