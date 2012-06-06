##
## JRuby adaptation layer for Kyoto Cabinet Java interface.
##
## This requires the kyotocabinet-java library to be installed. This
## library builds a jar file, which should be on the Ruby $LOAD_PATH,
## the Java CLASSPATH, or in /usr/local/lib/ where it is placed by the
## default 'make install' target in kyotocabinet-java.
##
## kyotocabinet-java also builds a necessary JNI library,
## libjkyotocabinet.{jnilib,dylib}. This cannot be added to
## java.library.path at runtime, so either its installation path must
## be added to java.library.path with a JVM option like
## -Djava.library.path=/some/dir, or it must be copied to somewhere
## already on java.library.path.
##

unless RUBY_PLATFORM == 'java'
  raise "This library requires JRuby!"
end

require 'java'

begin
  Java::Kyotocabinet::DB
rescue NameError
  # need to require kyotocabinet.jar
  if File.exist? '/usr/local/lib/kyotocabinet.jar'
    require '/usr/local/lib/kyotocabinet.jar'
    begin
      Java::Kyotocabinet::DB
    rescue NameError
      raise "Kyoto Cabinet classes could not be loaded, probably due to JNI library linking problems. Verify that libjkyotocabinet is in a directory on java.library.path (#{java.lang.System.getProperty('java.library.path')})."
    end
  else
    raise "kyotocabinet.jar could not be found!"
  end
end

module KyotoCabinet
  module Adaptation
    BYTE_ARRAY = [1].to_java(:byte).java_class

    def self.with_method_handle(mname, *args)
      yield self.java_method(mname, args)
    end
  end
end

module Java::Kyotocabinet
  class Cursor
    include KyotoCabinet::Adaptation

    with_method_handle(:get, BYTE_ARRAY, BYTE_ARRAY) do |m|
      def get(step=false)
        r = m.call(step)
        if r
          return [String.from_java_bytes(r[0]),
                  String.from_java_bytes(r[1])]
        else
          return nil
        end
      end
    end

    with_method_handle(:get_key, java.lang.Boolean) do |m|
      def get_key(step=false)
        r = m.call(step)
        if r
          return String.from_java_bytes(r)
        else
          return nil
        end
      end
    end

    with_method_handle(:jump, BYTE_ARRAY) do |m|
      def jump(key)
        m.call(key.to_java_bytes)
      end
    end
  end # class Cursor

  class DB
    include KyotoCabinet::Adaptation

    alias_method :_match_prefix, :match_prefix
    with_method_handle(:match_prefix, java.lang.String, java.lang.Integer) do |m|
      def match_prefix(prefix, limit=-1)
        m.call(prefix, limit)
      end
    end

    with_method_handle(:get, BYTE_ARRAY) do |m|
      def get(key)
        m.call(key)
      end
      alias_method :[], :get
    end

    with_method_handle(:set, BYTE_ARRAY, BYTE_ARRAY) do |m|
      def set(k, v)
        m.call(k.to_java_bytes, v.to_java_bytes)
      end
    end

    def cursor_process
      cur = self.cursor()
      begin
        yield cur
      ensure
        cur.disable()
      end
    end
  end # class DB
end

module KyotoCabinet
  import Java::Kyotocabinet
end
