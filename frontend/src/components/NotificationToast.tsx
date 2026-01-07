'use client';

import { useEffect } from 'react';
import { useStore, type Notification } from '@/lib/store';
import { X, CheckCircle, AlertCircle, Info, ArrowRight } from 'lucide-react';

const AUTO_DISMISS_MS = 10000;

function NotificationItem({ notification }: { notification: Notification }) {
  const { dismissNotification, setActiveTab } = useStore();

  useEffect(() => {
    const timer = setTimeout(() => {
      dismissNotification(notification.id);
    }, AUTO_DISMISS_MS);

    return () => clearTimeout(timer);
  }, [notification.id, dismissNotification]);

  const handleAction = () => {
    if (notification.action) {
      setActiveTab(notification.action.tab);
      dismissNotification(notification.id);
    }
  };

  const iconMap = {
    success: <CheckCircle className="w-5 h-5 text-green-400" />,
    error: <AlertCircle className="w-5 h-5 text-red-400" />,
    info: <Info className="w-5 h-5 text-blue-400" />,
  };

  const bgMap = {
    success: 'bg-green-900/90 border-green-700',
    error: 'bg-red-900/90 border-red-700',
    info: 'bg-blue-900/90 border-blue-700',
  };

  return (
    <div
      className={`${bgMap[notification.type]} border rounded-lg p-4 shadow-lg backdrop-blur-sm animate-slide-in`}
      role="alert"
    >
      <div className="flex items-start gap-3">
        {iconMap[notification.type]}
        <div className="flex-1 min-w-0">
          <p className="font-medium text-white">{notification.title}</p>
          <p className="text-sm text-gray-300 mt-0.5">{notification.message}</p>
          {notification.action && (
            <button
              onClick={handleAction}
              className="mt-2 inline-flex items-center gap-1.5 px-3 py-1.5 text-sm font-medium bg-white/10 hover:bg-white/20 rounded transition text-white"
            >
              {notification.action.label}
              <ArrowRight className="w-4 h-4" />
            </button>
          )}
        </div>
        <button
          onClick={() => dismissNotification(notification.id)}
          className="text-gray-400 hover:text-white transition"
          aria-label="Dismiss notification"
        >
          <X className="w-5 h-5" />
        </button>
      </div>
    </div>
  );
}

export function NotificationToast() {
  const { notifications } = useStore();

  if (notifications.length === 0) return null;

  return (
    <div className="fixed bottom-4 right-4 z-50 flex flex-col gap-2 max-w-sm">
      {notifications.map((notification) => (
        <NotificationItem key={notification.id} notification={notification} />
      ))}
    </div>
  );
}
