

# plot = figure(plot_width=800*_mult,
#               plot_height=400*_mult,
#               # background_fill_color='black',
#               # border_fill_color='black',
#               y_range=(0, max_intensity),
#               tools=TOOLS)

# hover_quads = HoverTool(renderers = [groups],
#                         tooltips=[('intensity', "@top{0,0}"),
#                                   ('m/z', "[@left{0,0.000}, @right{0,0.000}]")],
#                         mode='vline')
# hover_invisible = HoverTool(renderers = [invisible_buffers],
#                             tooltips=[('intensity', "@top{0,0}"),
#                                       ('m/z', "@x{0,0.000}")],
#                             mode='vline')
# plot.add_tools(hover_invisible, hover_bars, hover_squares, hover_quads)

# def get_buffers(mz_L, mz_R, max_buffer_length=.45):
#     r_prev = 0  # guard: 0th peak is a phoney
#     prev_width = 0
#     for i, (l, r) in enumerate(zip(mz_L, mz_R)):
#         buffer_length = min((l- r_prev)/2, max_buffer_length)
#         if i is not 0:
#             yield r_prev - prev_width/2, r_prev + buffer_length # previous right buffer
#         width = r - l
#         yield l - buffer_length, l + width/2  # current left buffer
#         r_prev = r
#         prev_width = width
#     yield r_prev, r_prev + max_buffer_length  # last right buffer
# hover_bars = HoverTool(renderers = [experimental_bars],
#                        tooltips=[('intensity', "@top{0,0}"),
#                                  ('m/z', "[@left{0,0.000}, @right{0,0.000}]")])
#
# hover_bars_e = HoverTool(renderers = [estimated_bars],
#                           tooltips=[('intensity', "@top{0,0}"),
#                                     ('m/z', "[@left{0,0.000},@right{0,0.000}]")])
# def get_buffers(mz_L, mz_R, max_length=.45):
#     L = []
#     R = []
#     r_prev = 0  # guards: 0th peak is a phoney
#     for i, (l, r) in enumerate(zip(mz_L, mz_R)):
#         tol = min((l - r_prev)/2, max_length)
#         if i:
#             R.append(r_prev + tol)
#         L.append(l - tol)
#         r_prev = r
#     R.append(r + tol)
#     return L, R
# not symmetric
# tolerance = .03
# experimental_bars = plot.vbar(x=masstodon.spectrum.mz,
#                               top=masstodon.spectrum.intensity,
#                               width=tolerance,
#                               color='black',
#                               alpha=.2)
# experimental_squares = plot.square(x=masstodon.spectrum.mz,
#                                    y=masstodon.spectrum.intensity,
#                                    size=5,
#                                    color='black',
#                                    alpha=.2)

# hover_invisible = HoverTool(renderers=[invisible_buffers],
#                             tooltips=[('intensity', "@top{0,0.000}"),
#                                       ('m/z', "@left{0,0.000}")],
#                             mode='vline')

def get_buffer(r_prev, l, r, l_next, max_length=.5):
    """Create one buffer."""
    tol = min(max_length, (l-r_prev)/2, (l_next-r)/2)
    return l - tol, r + tol

def buffer_args(L, R, max_length=.5):
    yield -infinity, L[0], R[0], L[1], max_length
    for i in range(1,len(L)-1):
        yield R[i-1], L[i], R[i] , L[i+1], max_length
    yield R[-2], L[-1], R[-1], infinity, max_length

def get_buffers(L, R, max_length=.5):
    """Get left and right ends of buffers."""
    return list(zip(*(get_buffer(*args) for args in buffer_args(mz_L, mz_R))))

## purely functional
# buffers_L, buffers_R = list(zip(*map(lambda arg: get_buffer(*arg),
#                                      buffer_args(mz_L, mz_R))))

## Pythonic
# buffers_L, buffers_R = list(zip(*(get_buffer(*args) for args in
#                                   buffer_args(mz_L, mz_R))))
