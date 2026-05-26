import numpy as np
import torch
import torch.optim as optim

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from torch.utils.data import DataLoader

from model import BinaryClassifier
from dataset import MergerDataset
from evaluate import evaluate
from utils import set_seed

#Setup
set_seed()

device = torch.device(
    "cuda" if torch.cuda.is_available() else "cpu"
)

print("Using device:", device)

#Loading Data
X = np.load("data/conc_features.npy")
y = np.load("data/conc_labels.npy")

print("X shape:", X.shape)
print("y shape:", y.shape)

#Normalize features
scaler = StandardScaler()

X = scaler.fit_transform(X)

#Train/Test Split
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42
)

train_dataset = MergerDataset(
    X_train,
    y_train
)

test_dataset = MergerDataset(
    X_test,
    y_test
)

#Initializing Dataloaders
train_loader = DataLoader(
    train_dataset,
    batch_size=512,
    shuffle=True
)

test_loader = DataLoader(
    test_dataset,
    batch_size=512,
    shuffle=False
)

#Loading Model
input_size = X.shape[1]

model = BinaryClassifier(
    input_size=input_size
).to(device)

#Loss and Optimizer
criterion = torch.nn.BCEWithLogitsLoss()

optimizer = optim.Adam(
    model.parameters(),
    lr=0.001
)

#Training Loop
epochs = 300

for epoch in range(epochs):

    model.train()

    running_loss = 0.0

    for X_batch, y_batch in train_loader:

        X_batch = X_batch.to(device)
        y_batch = y_batch.to(device)

        optimizer.zero_grad()

        outputs = model(X_batch)

        loss = criterion(
            outputs,
            y_batch
        )

        loss.backward()

        optimizer.step()
        running_loss += loss.item()

    avg_loss = running_loss / len(train_loader)

    accuracy = evaluate(
        model,
        test_loader,
        device
    )

    print(
        f"Epoch {epoch+1}/{epochs} | "
        f"Loss: {avg_loss:.4f} | "
        f"Accuracy: {accuracy:.4f}"
    )

#Model Saved
torch.save(
    model.state_dict(),
    "model.pth"
)

print("Model saved.")
