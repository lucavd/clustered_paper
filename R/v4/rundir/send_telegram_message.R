send_telegram_message <- function(message) {
  bot <- telegram.bot::Bot(token = Sys.getenv("TELEGRAM_BOT_TOKEN"))
  chat_id <- Sys.getenv("TELEGRAM_CHAT_ID")
  bot$sendMessage(chat_id = chat_id, text = message)
}
