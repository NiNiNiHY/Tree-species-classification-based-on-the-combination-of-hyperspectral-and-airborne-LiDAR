import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
import shap
from matplotlib import pyplot as plt

# from sklearn.linear_model import RandomizedLasso


# plt.rcParams.update({'figure.figsize': (12.0, 8.0)})
# plt.rcParams.update({'font.size': 14})
#
# boston = load_boston()
# X = pd.DataFrame(boston.data, columns=boston.feature_names)
# y = boston.target

#data = pd.read_excel('E:\\LIDAR\\project\\FeatureOfTree\\X.xlsx')
data = pd.read_excel('E:\\LIDAR\\project\\test40\\X.xlsx')
X = pd.DataFrame(data)

#输入特征组合，进行筛选，全部分析则注释
#X = X[['1', '6', '71',  '160',  'hmean', 'R90', 'SH' ,'Imean','R50' ,'f12','CF10']] #法国数据
#X = X[['5', '14', '195',  '205',  '426','CF5','CF6', 'R10', 'Isum' ,'Imean','Hmean','f8']] #neon数据

print(X)

#数据：
#data1 = pd.read_excel('E:\\LIDAR\\project\\FeatureOfTree\\Y.xlsx')  #法国数据
data1 = pd.read_excel('E:\\LIDAR\\project\\test40\\Y.xlsx')  #neon数据

y = np.array(data1)
y = y.ravel()
print(y)


print(X.shape)
print(y.shape)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)
rf = RandomForestRegressor(n_estimators=1000)  # 创建随机森林
rf.fit(X_train, y_train)  # 进行训练

print(rf.score(X_test, y_test))

rf.feature_importances_

print(X)
print(type(X))
print(type(y))
# print(type(boston))
sorted_idx = rf.feature_importances_.argsort()  # 对特征按照重要性进行排序

b = [sorted_idx >= 0]
c = sorted_idx[tuple(b)]
plt.barh(X_train.columns[c], rf.feature_importances_[c])  # 画图 纵轴为 特征名称 值为重要性
plt.xlabel("Random Forest Feature Importance")
plt.show()
print(c)
print(sorted_idx)
print(type(sorted_idx))

