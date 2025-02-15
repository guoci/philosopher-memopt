// Package cmd Slack top level command
package cmd

import (
	"errors"
	"os"

	"philosopher/lib/msg"
	"philosopher/lib/sla"
	"philosopher/lib/sys"

	"github.com/spf13/cobra"
)

// if you are here that's probably because you are curious on how Philosopher deals with the Slack user data.
// If you inspect the following code you will see that the organization is quite different from the other classes,
// all the information collected from the command-line execution is NOT internalized, and it's not added to the Meta
// data structure, which gathers the users input, so your private data is not stored and shared!

var name string
var direct string
var token string
var message string
var channel string

//var log string

// slackCmd represents the slack command
var slackCmd = &cobra.Command{
	Use:   "slack",
	Short: "Push notifications to Slack",
	Run: func(cmd *cobra.Command, args []string) {

		if len(token) < 1 {
			msg.InputNotFound(errors.New("you need to specify your token in order to push a notification"), "error")
		}

		sla.Run(name, direct, token, message, channel)
	},
}

func init() {
	if len(os.Args) > 1 && os.Args[1] == "slack" {

		m.Restore(sys.Meta())

		slackCmd.Flags().StringVarP(&name, "username", "", "Philosopher", "specify the name of the bot")
		slackCmd.Flags().StringVarP(&direct, "direct", "", "", "send a direct message to a user ID")
		slackCmd.Flags().StringVarP(&token, "token", "", "", "specify the Slack API token")
		slackCmd.Flags().StringVarP(&message, "message", "", "", "specify the text of the message to send")
		slackCmd.Flags().StringVarP(&channel, "channel", "", "", "specify the channel name")
		//slackCmd.Flags().StringVarP(&log, "log", "", "", "upload a log file")
	}

	RootCmd.AddCommand(slackCmd)
}
