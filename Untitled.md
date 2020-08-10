```python
#!/usr/bin/env python
# coding: utf-8

from sklearn import svm
import numpy as np
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import metrics
import matplotlib.pyplot as plt
```


```python

# ====================================
# STEP 1: read the training and testing data.
# Do not change any code of this step.

# specify path to training data and testing data
data_x = pd.read_csv('X_cancer.csv')
label_y = pd.read_csv("y_cancer.csv")

```


```python
scaler = StandardScaler().fit(data_x)
scaled_x = scaler.transform(data_x)
```


```python
type(scaled_x)
```




    numpy.ndarray




```python
label_y = label_y.to_numpy()
```


```python
type(label_y)
```




    numpy.ndarray




```python
X_train, X_test, y_train, y_test = train_test_split(
... scaled_x, label_y, test_size=0.33, random_state=42)
```


```python
#y_train
```


```python

# ====================================
# STEP 3: train model.
# Please modify the code in this step.
print("---train")
model = svm.SVC(kernel='poly',gamma=1.0) # this line should be changed
model.fit(X_train, y_train)

```

    ---train


    /opt/anaconda3/lib/python3.7/site-packages/sklearn/utils/validation.py:724: DataConversionWarning: A column-vector y was passed when a 1d array was expected. Please change the shape of y to (n_samples, ), for example using ravel().
      y = column_or_1d(y, warn=True)





    SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0,
        decision_function_shape='ovr', degree=3, gamma=1.0, kernel='poly',
        max_iter=-1, probability=False, random_state=None, shrinking=True,
        tol=0.001, verbose=False)




```python

# ====================================
# STEP3: evaluate model
# Don't modify the code below.
print("---evaluate")
print(" number of support vectors: ", model.n_support_)
acc = model.score(X_test, y_test)
print("acc:", acc)

```

    ---evaluate
     number of support vectors:  [34 22]
    acc: 0.9468085106382979



```python
y_pred_class = model.predict(X_test)
```


```python
#y_pred_class
```


```python
#y_test
```


```python
nnn=np.count_nonzero(np.absolute(y_test-y_pred_class))
```


```python
ddd=y_test.shape[0]
```


```python
ddd
```




    188




```python
nnn
```




    16616




```python
confusion = metrics.confusion_matrix(y_test, y_pred_class)
print(confusion)
```

    [[112   4]
     [  6  66]]



```python
TP = confusion[1, 1]
TN = confusion[0, 0]
FP = confusion[0, 1]
FN = confusion[1, 0]
```


```python
FP
```




    4




```python
# use float to perform true division, not integer division
print((TP + TN) / float(TP + TN + FP + FN))
print(metrics.accuracy_score(y_test, y_pred_class))
```

    0.9468085106382979
    0.9468085106382979



```python
# IMPORTANT: first argument is true values, second argument is predicted probabilities

# we pass y_test and y_pred_prob
# we do not use y_pred_class, because it will give incorrect results without generating an error
# roc_curve returns 3 objects fpr, tpr, thresholds
# fpr: false positive rate
# tpr: true positive rate
fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred_class)

plt.plot(fpr, tpr)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.rcParams['font.size'] = 12
plt.title('ROC curve for diabetes classifier')
plt.xlabel('False Positive Rate (1 - Specificity)')
plt.ylabel('True Positive Rate (Sensitivity)')
plt.grid(True)
```


![png](output_21_0.png)



```python

```
