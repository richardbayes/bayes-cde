function lgpdens_dashed(x, varargin) 
  ip = inputParser;
  ip.addRequired('x');
  ip.addParameter('bounded', [0 0]);
  ip.addParameter('percentiles', [2.5, 97.5]);
  ip.parse(x, varargin{:});
  x = ip.Results.x;
  bounded = ip.Results.bounded;
  percentiles = ip.Results.percentiles;
  [p, pq, xt] = lgpdens(x, 'bounded', bounded, 'percentiles', percentiles);
  newplot
  line(xt, p, 'linewidth', 2, 'Color', 'black', 'LineStyle', '--');
  line(xt, pq(:, 1), 'Color', 'black', 'LineStyle', '--');
  line(xt, pq(:, 2), 'Color', 'black', 'LineStyle', '--');
end