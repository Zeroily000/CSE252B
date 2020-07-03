import numpy as np
from PIL import Image, ImageOps
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import argparse
from scipy import signal, ndimage, interpolate
from typing import List

def featureDetector(img: np.ndarray, win_det: int, win_nms:int, feat_th: float) -> List[np.ndarray]:
    gray = np.array(ImageOps.grayscale(img))/255
    height, width = gray.shape

    # Compute gradients
    kernel_derivative = np.array([-1, 8, 0, -8, 1]) / 12
    Ix = ndimage.convolve1d(input=gray, weights=kernel_derivative, axis=1, mode='reflect')
    Iy = ndimage.convolve1d(input=gray, weights=kernel_derivative, axis=0, mode='reflect')

    # Compute the smaller eigenvalue of the gradient matrix
    kernel_window = np.ones(shape=(win_det, win_det)) / win_det ** 2
    IxIx = Ix * Ix
    IyIy = Iy * Iy
    IxIy = Ix * Iy
    IxIx_win_sum = signal.convolve2d(in1=IxIx, in2=kernel_window, mode='same', boundary='fill', fillvalue=0)
    IyIy_win_sum = signal.convolve2d(in1=IyIy, in2=kernel_window, mode='same', boundary='fill', fillvalue=0)
    IxIy_win_sum = signal.convolve2d(in1=IxIy, in2=kernel_window, mode='same', boundary='fill', fillvalue=0)

    grad_mat = np.zeros((height, width, 2, 2))
    grad_mat[:, :, 0, 0] = IxIx_win_sum
    grad_mat[:, :, 0, 1] = IxIy_win_sum
    grad_mat[:, :, 1, 0] = IxIy_win_sum
    grad_mat[:, :, 1, 1] = IyIy_win_sum
    eigen_value_min = np.linalg.eigvals(grad_mat)[:, :, -1]

    # Thresholding
    eigen_value_min[eigen_value_min < feat_th] = 0

    # Non-maximum suppression
    corners = np.logical_and(eigen_value_min != 0, eigen_value_min == ndimage.maximum_filter(eigen_value_min, (win_nms, win_nms)))

    # Subpixels
    x, y = np.meshgrid(np.arange(width), np.arange(height))

    xIxIx_win_sum = signal.convolve2d(in1=x*IxIx, in2=kernel_window, mode='same', boundary='fill', fillvalue=0)
    yIyIy_win_sum = signal.convolve2d(in1=y*IyIy, in2=kernel_window, mode='same', boundary='fill', fillvalue=0)
    xIxIy_win_sum = signal.convolve2d(in1=x*IxIy, in2=kernel_window, mode='same', boundary='fill', fillvalue=0)
    yIxIy_win_sum = signal.convolve2d(in1=y*IxIy, in2=kernel_window, mode='same', boundary='fill', fillvalue=0)

    keypoints = []
    for i in range(height):
        for j in range(width):
            if corners[i, j]:
                A = np.array([[IxIx_win_sum[i, j], IxIy_win_sum[i, j]],
                              [IxIy_win_sum[i, j], IyIy_win_sum[i, j]]])
                b = np.array([xIxIx_win_sum[i, j] + yIxIy_win_sum[i, j],
                              xIxIy_win_sum[i, j] + yIyIy_win_sum[i, j]])
                keypoints.append(np.linalg.inv(A).dot(b))
    return keypoints


def featureDescriptor(img: np.ndarray, keypoints: List[np.ndarray], win: int) -> List[np.ndarray]:
    gray = np.array(ImageOps.grayscale(img))/255
    height, width = gray.shape
    f = interpolate.interp2d(x=np.arange(width), y=np.arange(height), z=gray, kind='linear')

    descriptors = []
    for i, feat in enumerate(keypoints):
        x = np.linspace(feat[0] - win//2, feat[0] + win//2, win)
        y = np.linspace(feat[1] - win//2, feat[1] + win//2, win)
        descriptors.append(f(x, y).flatten())
    return descriptors


def featureMatcher(descl: List[np.ndarray], descr: List[np.ndarray], simi_th: float, dist_th: float):
    matched_points = []
    similarity_scores = corrCoef(np.array(descl), np.array(descr))
    while np.max(similarity_scores) > simi_th:
        i, j = np.unravel_index(np.argmax(similarity_scores), similarity_scores.shape)
        best_match_score = similarity_scores[i, j]
        similarity_scores[i, j] = -1
        next_best_match_scores = max(similarity_scores[i, :].max(), similarity_scores[:, j].max())

        if (1 - best_match_score) < (1 - next_best_match_scores)*dist_th:
            matched_points.append([i, j])
        similarity_scores[i, :] = -1
        similarity_scores[:, j] = -1
    return matched_points


def corrCoef(A: np.ndarray, B: np.ndarray):
    A_norm = A - A.mean(axis=1, keepdims=True)
    B_norm = B - B.mean(axis=1, keepdims=True)

    return A_norm.dot(B_norm.T) / np.sqrt(np.sum(A_norm*A_norm, axis=1, keepdims=True).dot(np.sum(B_norm*B_norm, axis=1, keepdims=True).T))


def plotFeatures(ax, keypoints, win):
    ax.scatter(np.array(keypoints).T[0], np.array(keypoints).T[1], s=1, c='r')
    for keypoint in keypoints:
        ax.add_patch(Rectangle(keypoint - win/2, win, win, fill=False))


def plotMatches(ax0, ax1, keypoints0, keypoints1, matches):
    for i, j in matches:
        ax0.plot([keypoints0[i][0], keypoints1[j][0]], [keypoints0[i][1], keypoints1[j][1]], '-r')
        ax1.plot([keypoints0[i][0], keypoints1[j][0]], [keypoints0[i][1], keypoints1[j][1]], '-r')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Corner detection')
    parser.add_argument('--left-image', '-l',
                        dest='imgl',
                        type=str,
                        required=True)
    parser.add_argument('--right-image', '-r',
                        dest='imgr',
                        type=str,
                        required=True)
    parser.add_argument('--win-det',
                        dest='win_det',
                        type=int,
                        required=False,
                        default=11)
    parser.add_argument('--win-nms',
                        dest='win_nms',
                        type=int,
                        required=False,
                        default=7)
    parser.add_argument('--feat-th',
                        dest='feat_th',
                        type=float,
                        required=False,
                        default=0.0009)
    parser.add_argument('--simi-th',
                        dest='simi_th',
                        type=float,
                        required=False,
                        default=0.5)
    parser.add_argument('--dist-th',
                        dest='dist_th',
                        type=float,
                        required=False,
                        default=0.6)
    args = parser.parse_args()

    # Load images
    imgl = Image.open(fp=args.imgl)
    imgr = Image.open(fp=args.imgr)

    # Detect features
    kptl = featureDetector(img=imgl, win_det=args.win_det, win_nms=args.win_nms, feat_th=args.feat_th)
    kptr = featureDetector(img=imgr, win_det=args.win_det, win_nms=args.win_nms, feat_th=args.feat_th)

    # Compute descriptors
    descl = featureDescriptor(img=imgl, win=args.win_det, keypoints=kptl)
    descr = featureDescriptor(img=imgr, win=args.win_det, keypoints=kptr)

    # Feature matching
    matched_features = featureMatcher(descl=descl, descr=descr, simi_th=args.simi_th, dist_th=args.dist_th)
    print(f'Matched {len(matched_features)} features')

    # Plot
    fig, axes = plt.subplots(nrows=1, ncols=2)
    axes[0].imshow(imgl)
    axes[1].imshow(imgr)
    plotFeatures(ax=axes[0], keypoints=kptl, win=args.win_det)
    plotFeatures(ax=axes[1], keypoints=kptr, win=args.win_det)
    plotMatches(ax0=axes[0], ax1=axes[1], keypoints0=kptl, keypoints1=kptr, matches=matched_features)
    axes[0].set_title(f'Detected {len(kptl)} features')
    axes[1].set_title(f'Detected {len(kptr)} features')
    fig.show()


