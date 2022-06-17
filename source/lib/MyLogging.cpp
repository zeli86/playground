

#include "MyLogging.h"

namespace logging = boost::log;
namespace attrs = boost::log::attributes;
namespace expr = boost::log::expressions;
namespace keywords = boost::log::keywords;
namespace sinks = boost::log::sinks;


void initLogging(const std::string& sLogFilename)
{
  typedef sinks::synchronous_sink<sinks::text_ostream_backend> text_sink;

  {
    boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();
    sink->locked_backend()->add_stream(boost::shared_ptr<std::ostream>(&std::cout));

// flush
    sink->locked_backend()->auto_flush(true);

    // format sink
    // sink->set_formatter
    // (
    //   /// TODO add your desired formatting
    // );

    // filter
    // TODO add any filters

    // register sink
    logging::core::get()->add_sink(sink);
  }

  {
// create sink to logfile
    boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();
    sink->locked_backend()->add_stream(boost::make_shared<std::ofstream>(sLogFilename.c_str()));

    // flush
    sink->locked_backend()->auto_flush(true);

    // format sink
    // sink->set_formatter
    // (
    //     /// TODO add your desired formatting
    // );

    // filter
    sink->set_filter(logging::trivial::severity >= logging::trivial::info);

    // register sink
    logging::core::get()->add_sink(sink);
  }

  // Also let's add some commonly used attributes, like timestamp and record counter.
  logging::add_common_attributes();
  logging::core::get()->add_thread_attribute("Scope", attrs::named_scope());
}

void disableLogging()
{
  logging::core::get()->set_logging_enabled(false);
}