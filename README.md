# ggimg
Poor man's 2d and 3d image operations

- Header-only
- No 3rd party libraries

Below are some of the available operations. Check the header file for the complete set of operations.

---

* original

<a href="https://i.imgur.com/2iz4hZD.png" target="_blank">![original](https://i.imgur.com/2iz4hZD.png)</a>

* ggimg::rgb_to_luma601_2d()

<a href="https://i.imgur.com/TG44h6G.png" target="_blank">![rgb_to_luma601_2d](https://i.imgur.com/TG44h6G.png)</a>

* ggimg::rgb_to_luma709_2d()

<a href="https://i.imgur.com/PwKGOlN.png" target="_blank">![rgb_to_luma709_2d](https://i.imgur.com/PwKGOlN.png)</a>

* ggimg::rgb_to_gray_2d()

<a href="https://i.imgur.com/YpYUQKD.png" target="_blank">![rgb_to_gray_2d](https://i.imgur.com/YpYUQKD.png)</a>

* ggimg::normalize_2d(dmin = 100, dmax = 200)

<a href="https://i.imgur.com/hUeKw8k.png" target="_blank">![normalize_2d](https://i.imgur.com/hUeKw8k.png)</a>

* ggimg::normalize_hist_2d()

<a href="https://i.imgur.com/SID2yet.png" target="_blank">![normalize_hist_2d](https://i.imgur.com/SID2yet.png)</a>

* ggimg::gradient_sobel_2d(mode = 1)

<a href="https://i.imgur.com/VeyG5N0.png" target="_blank">![gradient_sobel_2d](https://i.imgur.com/VeyG5N0.png)</a>

* ggimg::gradient_sobel_2d(mode = 2)

<a href="https://i.imgur.com/6cPw22Z.png" target="_blank">![gradient_sobel_2d](https://i.imgur.com/6cPw22Z.png)</a>

* ggimg::gradient_sobel_2d(mode = 0)

<a href="https://i.imgur.com/1cw11Fv.png" target="_blank">![gradient_sobel_2d](https://i.imgur.com/1cw11Fv.png)</a>

* ggimg::gaussian_filter_2d(sigma = 3.0f)

<a href="https://i.imgur.com/GofjQqU.png" target="_blank">![gaussian_filter_2d](https://i.imgur.com/GofjQqU.png)</a>

* ggimg::median_filter_2d(k = 5)

<a href="https://i.imgur.com/JxftfXY.png" target="_blank">![median_filter_2d](https://i.imgur.com/JxftfXY.png)</a>

* ggimg::scale_nn_2d(sx = 0.33f, sy = 0.66f)

<a href="https://i.imgur.com/RyVTdak.png" target="_blank">![scale_nn_2d](https://i.imgur.com/RyVTdak.png)</a>

* ggimg::scale_li_2d(sx = 0.33f, sy = 0.66f)

<a href="https://i.imgur.com/99qGoZ5.png" target="_blank">![scale_li_2d](https://i.imgur.com/99qGoZ5.png)</a>

