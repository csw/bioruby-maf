When /^slice the resulting blocks according to the given interval$/ do
  # @blocks and @interval
  @blocks = @blocks.collect { |b| b.slice(@interval) }
end
