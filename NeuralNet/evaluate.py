import torch

def evaluate(model, loader, device):

    model.eval()

    correct = 0
    total = 0

    with torch.no_grad():

        for X_batch, y_batch in loader:

            X_batch = X_batch.to(device)
            y_batch = y_batch.to(device)

            outputs = model(X_batch)

            probs = torch.sigmoid(outputs)

            preds = (probs > 0.5).float()

            correct += (preds == y_batch).sum().item()

            total += y_batch.size(0)

    accuracy = correct / total

    return accuracy
