import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def show_images_grid(image_paths, dpi=300):
    fig, axs = plt.subplots(2, 3, figsize=(15, 10), dpi=dpi)  # 創建一個3x3的圖表，每個單元可以放一張圖片
    axs = axs.ravel()  # 將二維的axes陣列變為一維，便於迭代
    # fig.suptitle('With vs. Without Filter CCCHCCHH', fontsize=20, y=0.95)
    for i, img_path in enumerate(image_paths):
        img = mpimg.imread(img_path)  # 讀取圖片
        axs[i].imshow(img)  # 將圖片顯示在對應的子圖上
        # axs[i].set_title(f'Image {i+1}')  # 為每個子圖添加標題
        axs[i].axis('off')  # 不顯示坐標軸

    plt.tight_layout()  # 調整子圖間距
    plt.savefig('without_filter_CCCHCCHH.png')  # 將圖表儲存為圖片檔案

# 圖片檔案名稱列表
image_files = [f'./i{n}.png' for n in range(1, 7)]

# 呼叫函數顯示圖片
show_images_grid(image_files, dpi=500)